"""FASTQ-only pipeline: no Streamlit, clean library, lineage or original payload.

Deterministic progressive profile MSA (not POA). Gap/base uncertainty is retained
with alternatives for bounded candidate search. Public settings are explicit.
"""
from dataclasses import dataclass, asdict
from collections import Counter, defaultdict
from functools import lru_cache
import gzip
import io
import math
import hashlib
import numpy as np
import hedges_bridge as hb
import legacy_codec as lc
import framing
from process_model import rc


@dataclass(frozen=True)
class Settings:
    codec: str = "TVH-V5.2 HEDges"
    fbr: str = hb.FBR
    rbr: str = hb.RBR
    sid_mode: str = hb.SID_MODE_HAMMING
    primer_edits: int = 2
    si_edits: int = 1
    max_edits: int = 1
    min_q: int = 0
    min_reads: int = 1
    max_reads: int = 500
    support: float = 0.6
    method: str = "Quality-weighted voting"
    outlier_fraction: float = 0.25
    min_overlap: int = 25
    overlap_error: float = 0.15
    expected_strands: int = 0
    stored_bytes: int = 0
    dna_length: int = 0
    joint_si_body: bool = False
    rescue_uncertain: bool = True
    max_uncertain: int = 6
    max_candidates: int = 256
    search_seconds: float = 30.0


def validate(c):
    if not 0<=c.max_uncertain<=12 or not 1<=c.max_candidates<=4096 or not .1<=c.search_seconds<=600:
        raise ValueError('Invalid candidate search limits')
    if c.codec not in ("TVH-V5.2 HEDges","Simple Mapping","R0_B9","R1_B12","R2_B15","RINF_B16","RINF_SHUFFLED_B16"):
        raise ValueError("Unknown DNA encoding")
    if c.fbr != hb.FBR or c.rbr != hb.RBR:
        raise ValueError("This release uses the original app's primer pair")
    if c.sid_mode not in hb.SID_MODES or not 0<=c.primer_edits<=3 or not 0<=c.si_edits<=2 or not 0<=c.max_edits<=4:
        raise ValueError("Invalid primer/SI/decoder settings")
    if not .5<c.support<=1 or not 1<=c.min_reads<=c.max_reads<=500 or not 0<=c.min_q<=60:
        raise ValueError("Invalid consensus settings")
    if not 0<=c.expected_strands<=65536 or not 0<=c.stored_bytes<=100_000_000 or not 0<=c.dna_length<=10_000_000:
        raise ValueError("Invalid public framing lengths")
    return c


def load_settings(raw):
    obj=__import__('json').loads(raw)
    if not isinstance(obj,dict) or set(obj)-set(Settings.__dataclass_fields__):
        raise ValueError("Settings contain unsupported fields. No payload/reference/ground-truth fields are allowed.")
    return validate(Settings(**obj))


def parse(raw):
    if raw[:2]==b'\x1f\x8b':
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f: raw=f.read(100_000_001)
    if len(raw)>100_000_000: raise ValueError("FASTQ exceeds the 100 MB uncompressed limit")
    lines=raw.decode('ascii').splitlines()
    if len(lines)%4: raise ValueError("FASTQ must contain four lines per record")
    records=[]
    for i in range(0,len(lines),4):
        h,s,p,q=lines[i:i+4]; s=s.upper()
        if not h.startswith('@') or not p.startswith('+') or not s or len(s)!=len(q) or set(s)-set('ACGTN') or any(not 33<=ord(x)<=126 for x in q):
            raise ValueError(f"Invalid FASTQ record {i//4+1}")
        if len(s)>600: raise ValueError("Reads longer than 600 nt are outside this oligo pipeline")
        records.append({'name':h[1:],'seq':s,'q':[ord(x)-33 for x in q]})
    return records


@lru_cache(maxsize=16384)
def distance(a,b):
    if a==b:return 0
    row=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        new=[i]
        for j,y in enumerate(b,1): new.append(min(new[-1]+1,row[j]+1,row[j-1]+(x!=y)))
        row=new
    return row[-1]


@lru_cache(maxsize=4096)
def primer_hits(seq,primer,k):
    # Semi-global edit alignment. Free read flanks; primer must align in full.
    m,n=len(primer),len(seq)
    d=[[0]*(n+1) for _ in range(m+1)]
    for i in range(1,m+1): d[i][0]=i
    for i in range(1,m+1):
        for j in range(1,n+1): d[i][j]=min(d[i-1][j-1]+(primer[i-1]!=seq[j-1]),d[i-1][j]+1,d[i][j-1]+1)
    best=min(d[m])
    if best>k: return ()
    hits=set()
    for end in range(1,n+1):
        if d[m][end]!=best: continue
        # Enumerate tied starts, bounded by primer search window.
        pending=[(m,end)]; visited=set()
        while pending:
            i,j=pending.pop()
            if (i,j) in visited: continue
            visited.add((i,j))
            if i==0:
                if end-j>=m-k: hits.add((j,end,best))
                continue
            if j and d[i][j]==d[i-1][j-1]+(primer[i-1]!=seq[j-1]): pending.append((i-1,j-1))
            if d[i][j]==d[i-1][j]+1: pending.append((i-1,j))
            if j and d[i][j]==d[i][j-1]+1: pending.append((i,j-1))
    return tuple(sorted(hits))


def orient_trim(rec,c):
    hits=[]
    for seq,q in ((rec['seq'],rec['q']),(rc(rec['seq']),rec['q'][::-1])):
        forward=primer_hits(seq[:70],c.fbr,c.primer_edits)
        shift=max(0,len(seq)-100)
        reverse=primer_hits(seq[shift:],c.rbr,c.primer_edits)
        for fs,fe,fc in forward:
            for rs,re,rr in reverse:
                rs+=shift
                if rs>fe+5:
                    hits.append((fc+rr,seq[fe:rs],tuple(q[fe:rs])))
    if not hits: return None,'PRIMER_NOT_FOUND'
    score=min(h[0] for h in hits)
    best={(s,q) for cost,s,q in hits if cost==score}
    if len(best)!=1: return None,'PRIMER_AMBIGUOUS'
    seq,q=next(iter(best))
    return {'seq':seq,'q':list(q)},'OK'


def progressive_msa(reads):
    """Profile MSA with per-read base qualities; gap columns are retained.
    Read ordering is canonical and independent of FASTQ names/order.
    """
    observed=Counter(r['seq'] for r in reads)
    if len(observed)==1:
        return [[(r['seq'][i],r['q'][i]) for r in reads] for i in range(len(reads[0]['seq']))]
    guide=min(sorted(observed)[:12],key=lambda s:(sum(distance(s,t)*n for t,n in observed.items()),s))
    ordered=sorted(reads,key=lambda r:(distance(guide,r['seq']),-sum(r['q'])/len(r['q']),r['seq'],tuple(r['q'])))
    columns=[[(b,q)] for b,q in zip(ordered[0]['seq'],ordered[0]['q'])]
    used=1
    for r in ordered[1:]:
        seq,quals=r['seq'],r['q']; m,n=len(columns),len(seq)
        d=np.zeros((m+1,n+1),dtype=float); bt=np.zeros((m+1,n+1),dtype=np.int8)
        d[:,0]=-3*np.arange(m+1); d[0,:]=-3*np.arange(n+1)
        for i in range(1,m+1):
            counts=Counter(x for x,_ in columns[i-1]); denom=len(columns[i-1])
            non_gap=denom-counts['-']
            for j in range(1,n+1):
                match=(2*counts[seq[j-1]]-2*(non_gap-counts[seq[j-1]]-counts['N'])-3*counts['-'])/denom if seq[j-1]!='N' else -1
                values=(d[i-1,j-1]+match,d[i-1,j]-3*non_gap/denom,d[i,j-1]-3)
                ix=int(np.argmax(values)); d[i,j]=values[ix]; bt[i,j]=ix
        new=[]; i,j=m,n
        while i or j:
            t=int(bt[i,j]) if i and j else (1 if i else 2)
            if t==0: new.append(columns[i-1]+[(seq[j-1],quals[j-1])]); i-=1;j-=1
            elif t==1: new.append(columns[i-1]+[('-',0)]);i-=1
            else: new.append([('-',0)]*used+[(seq[j-1],quals[j-1])]);j-=1
        columns=new[::-1];used+=1
    return columns


def consensus(reads,c):
    ordered=sorted(reads,key=lambda r:(-sum(r['q'])/len(r['q']),r['seq'],tuple(r['q'])))[:c.max_reads]
    if not ordered: return '',[],{'used':0,'outliers':0,'uncertain':0}
    # Medoid selected from observed sequences only; all ties resolved canonically.
    unique=sorted({r['seq'] for r in ordered})
    candidates=unique[:12]
    counts=Counter(r['seq'] for r in ordered)
    medoid=min(candidates,key=lambda s:(sum(distance(s,t)*n for t,n in counts.items()),s))
    kept=[r for r in ordered if distance(medoid,r['seq'])/max(len(medoid),len(r['seq']),1)<=c.outlier_fraction]
    info={'used':len(kept),'outliers':len(ordered)-len(kept),'capped':max(0,len(reads)-len(ordered)),'uncertain':0}
    if len(kept)<c.min_reads: return '',[],info
    columns=progressive_msa(kept); out=[]; pos=[]
    for index,col in enumerate(columns):
        gap=sum(b=='-' for b,q in col)/len(col)
        weights=Counter()
        for b,q in col:
            if b not in 'ACGT': continue
            weights[b]+=1 if c.method=='Majority voting' else 1-10**(-min(q,60)/10)
        total=sum(weights.values()); best=max(weights.values(),default=0)
        winners=[b for b,v in weights.items() if abs(v-best)<1e-10]
        base_support=best/total if total else 0
        # Structural gap uncertainty cannot be hidden by choosing a nucleotide.
        if gap>=c.support: base='-'
        elif 1-gap>=c.support and len(winners)==1 and base_support>=c.support: base=winners[0]
        else: base='N';info['uncertain']+=1
        if base!='-': out.append(base)
        alternatives=[b for b in 'ACGT' if weights[b]>0]
        if not alternatives:alternatives=list('ACGT')
        if gap>0:alternatives.append('')
        pos.append({'column':index+1,'base':base,'base_support':base_support,'gap_fraction':gap,'reads':len(col),
                    'alternatives':alternatives,'base_weights':dict(weights)})
    return ''.join(out),pos,info


def merge_pair(a,b,c):
    """Indel-aware suffix-prefix overlap; quality disagreement becomes N.
    Standard paired-end orientation: R2 reverse-complemented before assembly.
    """
    x,y=a['seq'],rc(b['seq']); qx,qy=a['q'],b['q'][::-1]
    m,n=len(x),len(y)
    d=np.zeros((m+1,n+1)); bt=np.zeros((m+1,n+1),np.int8)
    d[0,:]=-3*np.arange(n+1)  # free prefix of R1, mandatory prefix of R2
    for i in range(1,m+1):
        for j in range(1,n+1):
            vals=(d[i-1,j-1]+(2 if x[i-1]==y[j-1] else -3),d[i-1,j]-3,d[i,j-1]-3)
            t=int(np.argmax(vals));d[i,j]=vals[t];bt[i,j]=t
    ends=[j for j in range(c.min_overlap,n+1)]
    if not ends:return None
    score=max(d[m,j] for j in ends)
    endpoints=[j for j in ends if d[m,j]==score]
    candidates=[]
    for end in endpoints:
        i,j=m,end; columns=[]; edits=0; matched=0
        while j:
            if not i: break
            t=int(bt[i,j])
            if t==0:
                columns.append(((x[i-1],qx[i-1]),(y[j-1],qy[j-1])));edits+=x[i-1]!=y[j-1];matched+=1;i-=1;j-=1
            elif t==1:columns.append(((x[i-1],qx[i-1]),('-',0)));edits+=1;i-=1
            else:columns.append((('-',0),(y[j-1],qy[j-1])));edits+=1;j-=1
        if j or matched<c.min_overlap or edits/max(1,len(columns))>c.overlap_error:continue
        s=list(x[:i]);q=list(qx[:i])
        for (v,qv),(w,qw) in reversed(columns):
            if v==w:s.append(v);q.append(min(qv,qw))
            elif '-' in (v,w):s.append('N');q.append(0)
            elif abs(qv-qw)>=6:s.append(v if qv>qw else w);q.append(max(qv,qw)-3)
            else:s.append('N');q.append(0)
        s.extend(y[end:]);q.extend(qy[end:]);candidates.append((''.join(s),tuple(q)))
    unique=set(candidates)
    if len(unique)!=1:return None
    seq,q=next(iter(unique));return {'seq':seq,'q':list(q)}


@lru_cache(maxsize=8192)
def si_candidates(prefix,mode,k):
    cb=hb.codebook(hb.strand_capacity(mode),mode)
    exact=[sid for sid,si in cb.items() if si==prefix[:8]]
    if exact:return tuple((sid,8) for sid in exact)
    if mode==hb.SID_MODE_NO_HAMMING:return ()
    hits=[]
    for length in range(max(1,8-k),min(len(prefix),8+k)+1):
        observed=prefix[:length]
        if 'N' in observed:continue
        for sid,si in cb.items():
            d=distance(observed,si)
            if d<=k:hits.append((d,sid,length))
    if not hits:return ()
    best=min(h[0] for h in hits)
    return tuple(sorted({(sid,length) for d,sid,length in hits if d==best}))


def assign(interior,c):
    seq,q=interior['seq'],interior['q']
    if c.codec!='TVH-V5.2 HEDges':
        if len(seq)<8 or set(seq[:8])-set('ACGT'):return None,'SI_UNREADABLE'
        sid=int(''.join(f'{"ACGT".index(b):02b}' for b in seq[:8]),2)
        hits=[(sid,8)]
    else:
        radius=c.si_edits if c.sid_mode==hb.SID_MODE_HAMMING else 0
        hits=list(si_candidates(seq[:8+radius],c.sid_mode,radius))
    if c.expected_strands: hits=[h for h in hits if h[0]<c.expected_strands]
    if not hits:return None,'SI_UNASSIGNED'
    if len(hits)>1 and c.joint_si_body and c.codec=='TVH-V5.2 HEDges':
        if len(hits)>8:return None,'SI_SEARCH_LIMIT'
        verified=[]
        for sid,length in hits:
            body=seq[length:]
            if 'N' in body or abs(len(body)-76)>c.max_edits:continue
            d=hb.decode_body(body,sid,c.max_edits,c.sid_mode)
            if d.status=='OK':verified.append((sid,length))
        hits=verified
    if len(hits)!=1:return None,'SI_AMBIGUOUS'
    sid,length=hits[0]
    return (sid,{'seq':seq[length:],'q':q[length:]}),'OK'


def run(r1,c=Settings(),r2=b'',progress=None,decode_file=True):
    validate(c);a=parse(r1);b=parse(r2) if r2 else None
    if not a:raise ValueError('FASTQ contains no reads')
    if b is not None:
        def name(s):
            s=s.split()[0];return s[:-2] if s.endswith(('/1','/2')) else s
        if len(a)!=len(b) or any(name(x['name'])!=name(y['name']) for x,y in zip(a,b)):
            raise ValueError('R1/R2 record counts or pair names do not match')
    groups=defaultdict(list);counts=Counter();audit=[];observed=[]
    for index,rec in enumerate(a):
        mates=[rec] if b is None else [rec,b[index]]
        valid=[r for r in mates if np.mean(r['q'])>=c.min_q]
        counts['quality_filtered']+=len(mates)-len(valid)
        calls=[]
        for r in valid:
            interior,status=orient_trim(r,c)
            if interior:
                hit,status=assign(interior,c)
                if hit:calls.append(hit)
            counts[status]+=1
        if len(calls)==2:
            if calls[0][0]!=calls[1][0]:counts['PAIR_SI_CONFLICT']+=1;calls=[]
            else:
                # One full-insert mate, chosen by observed quality only. Do not
                # count two correlated mates as independent cluster evidence.
                calls=[max(calls,key=lambda x:(sum(x[1]['q'])/len(x[1]['q']),x[1]['seq'],tuple(x[1]['q'])))]
                counts['PAIRED_FULL_INSERT_SELECTION']+=1
        elif not calls and len(valid)==2:
            merged=merge_pair(valid[0],valid[1],c)
            if merged:
                interior,status=orient_trim(merged,c)
                if interior:
                    hit,status=assign(interior,c)
                    if hit:calls=[hit];counts['PAIRS_MERGED']+=1
                counts['merged_'+status]+=1
            else:counts['PAIR_OVERLAP_FAILED']+=1
        if calls:
            sid,body=calls[0]
            if body['seq']:groups[sid].append(body)
            audit.append(sid)
            observed.append({'unit':index,'sid':sid,'sequence':body['seq']})
        else:audit.append(None)
        if progress and index%20==0:progress(.5*(index+1)/len(a))
    bodies={};all_bodies={};outcomes=[];positions=[]
    for index,sid in enumerate(sorted(groups)):
        seq,pos,info=consensus(groups[sid],c)
        all_bodies[sid]=seq
        if seq and 'N' not in seq:bodies[sid]=seq
        outcomes.append({'Strand ID':sid,'Reads':len(groups[sid]),'Consensus length':len(seq),
                         'Consensus status':'OK' if seq and 'N' not in seq else 'UNCERTAIN' if seq else 'INSUFFICIENT_READS',**info})
        positions.extend({'Strand ID':sid,**p} for p in pos)
        if progress:progress(.5+.5*(index+1)/len(groups))
    result={'status':'CONSENSUS_READY','data':None}
    result.update({'consensus':bodies,'all_consensus':all_bodies,'observed_bodies':observed,'consensus_rows':outcomes,'positions':positions,
                   'read_stats':{'input_units':len(a),'raw_reads':len(a)*(2 if b is not None else 1),
                                 'assigned_units':sum(v is not None for v in audit),**dict(counts)},
                   'settings':asdict(c),'assignments':audit})
    if decode_file:
        from candidate_search import decode
        result.update(decode(result,c))
    return result


def reconstruct(bodies,c,verified=None):
    expected=list(range(c.expected_strands)) if c.expected_strands else list(range(max(bodies,default=-1)+1))
    if not expected or sorted(bodies)!=expected:
        return {'status':'FAIL','error':'Missing or uncertain strand(s). No outer ECC is used.','data':None,'strand_results':[]}
    if c.codec=='TVH-V5.2 HEDges':
        decoded=[verified[sid] if verified is not None else hb.decode_body(bodies[sid],sid,c.max_edits,c.sid_mode) for sid in expected]
        if any(d.status!='OK' for d in decoded):
            return {'status':'FAIL','error':'HEDGES could not uniquely decode every strand.','data':None,'strand_results':[d.__dict__ for d in decoded]}
        bits=''.join(d.payload_bits for d in decoded)
        embedded_magic=''.join(f'{b:08b}' for b in framing.HEDGES_FILE_MAGIC)
        if bits.startswith(embedded_magic):
            status,chosen,candidates=framing.autonomous_reconstruct(bits)
            result={'status':status,'data':chosen.data if chosen else None,'format':chosen.format if chosen else None,
                    'extension':chosen.extension if chosen else '.bin','mime':chosen.mime if chosen else 'application/octet-stream',
                    'filename':chosen.filename if chosen else None,'domain':chosen.domain if chosen else None,
                    'padding_bits':chosen.padding_bits if chosen else None,
                    'strand_results':[d.__dict__ for d in decoded]}
        elif c.stored_bytes:
            end=c.stored_bytes*8; tail=bits[end:]
            if end>len(bits) or (tail and tail!='1'+'0'*(len(tail)-1)):
                return {'status':'FAIL','error':'Public file size does not match decoded framing','data':None,'strand_results':[d.__dict__ for d in decoded]}
            data=lc.bits_to_bytes(bits[:end]); fmt,ext,mime,strong,note=framing.detect_and_validate_exact(data)
            result={'status':'OK','data':data,'format':fmt,'extension':ext,'mime':mime,'strand_results':[d.__dict__ for d in decoded]}
        else:
            status,chosen,candidates=framing.autonomous_reconstruct(bits)
            result={'status':status,'data':chosen.data if chosen else None,'format':chosen.format if chosen else None,
                    'extension':chosen.extension if chosen else '.bin','mime':chosen.mime if chosen else 'application/octet-stream',
                    'strand_results':[d.__dict__ for d in decoded]}
    else:
        if not c.dna_length or not c.stored_bytes:
            return {'status':'FAIL','error':'For rule-based/Simple Mapping, provide public encoded DNA length and stored byte count.','data':None,'strand_results':[]}
        if any(len(bodies[i])!=76 for i in expected) or math.ceil(c.dna_length/76)!=len(expected):
            return {'status':'FAIL','error':'Legacy consensus lengths/layout do not match the public design. No silent padding or truncation repairs.','data':None,'strand_results':[]}
        dna=''.join(bodies[i] for i in expected)[:c.dna_length]
        try:data,_,_=lc.decode_dna(dna,c.codec,{'bits_len':8*c.stored_bytes,'bytes_len':c.stored_bytes})
        except Exception as exc:return {'status':'FAIL','error':str(exc),'data':None,'strand_results':[]}
        fmt,ext,mime,strong,note=framing.detect_and_validate_exact(data)
        result={'status':'OK','data':data,'format':fmt,'extension':ext,'mime':mime,'strand_results':[]}
    result['validation_note']='Decoded candidate; exact byte equality requires an external hash check. Unknown final missing strands cannot be excluded when expected count is unspecified.'
    if result.get('data') is not None:
        result['sha256']=hashlib.sha256(result['data']).hexdigest()
        restored=framing.restore_lossless_payload(result['data'])
        if restored:
            result.update(restored_data=restored['data'],restored_name=restored['name'],restored_mime=restored['mime'])
    return result
