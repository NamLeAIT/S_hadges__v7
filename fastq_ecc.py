from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip, io, math, random, re
from typing import Dict, List, Tuple, Any

BASES='ACGT'

def clean_dna(s:str)->str:
    return ''.join(c for c in str(s or '').upper() if c in BASES)


def parse_fastq(data: bytes|str)->List[Dict[str,Any]]:
    if isinstance(data,(bytes,bytearray)):
        raw=bytes(data)
        if raw[:2]==b'\x1f\x8b':
            raw=gzip.decompress(raw)
        text=raw.decode('utf-8-sig',errors='ignore')
    else:
        text=str(data)
    lines=[x.rstrip('\r') for x in text.splitlines() if x.strip()!='']
    if len(lines)%4: raise ValueError('FASTQ must contain complete four-line records')
    out=[]
    for i in range(0,len(lines),4):
        h,s,p,q=lines[i:i+4]
        if not h.startswith('@') or not p.startswith('+'):
            raise ValueError(f'Malformed FASTQ record {i//4+1}')
        s=clean_dna(s)
        if len(s)!=len(q): raise ValueError(f'FASTQ record {i//4+1}: sequence/quality length mismatch')
        out.append({'header':h[1:].strip(),'sequence':s,'quality':q,'phred':[max(0,ord(c)-33) for c in q]})
    return out


def sid_from_header(header:str)->int|None:
    pats=(r'(?:strand|sid|strand_id|strand-id)[_:=\- ]*(\d+)', r'(?:^|[|:/_\-])(\d+)(?:$|[|:/_\-])')
    for p in pats:
        m=re.search(p,str(header),flags=re.I)
        if m: return int(m.group(1))
    return None


def _qchar(q:int)->str:
    return chr(max(0,min(93,int(q)))+33)

Q_MIN = 2
Q_MAX = 40
STRUCTURAL_Q = 40


def _clamp_q(q: float) -> int:
    return max(Q_MIN, min(Q_MAX, int(round(float(q)))))


def _q_to_error_probability(q: float) -> float:
    return 10.0 ** (-float(q) / 10.0)


def _p_to_q(p: float) -> int:
    p = min(1.0, max(1e-12, float(p)))
    return _clamp_q(-10.0 * math.log10(p))


def make_qscores(
    length: int,
    rng: random.Random,
    profile: str = 'Fixed Q',
    q_mean: float = 30.0,
    q_std: float = 4.0,
    q_start: float = 35.0,
    q_end: float = 20.0,
    q_fixed: float = 30.0,
) -> List[int]:
    """App-ECC Phred profiles used by the Q-score-driven substitution model."""
    n=max(0,int(length)); profile=str(profile or 'Fixed Q')
    if n==0: return []
    if profile=='Normal Q distribution':
        return [_clamp_q(rng.gauss(float(q_mean), max(0.0,float(q_std)))) for _ in range(n)]
    if profile=='Linear quality decay':
        if n==1: return [_clamp_q(q_start)]
        return [_clamp_q(float(q_start)+(float(q_end)-float(q_start))*(i/(n-1))) for i in range(n)]
    return [_clamp_q(q_fixed)]*n


def _mutate_segment(
    seq:str,
    rng:random.Random,
    substitution_model:str,
    configured_substitution_rate:float,
    qscores:List[int],
    insertion_rate:float,
    deletion_rate:float,
)->Tuple[str,List[int],Dict[str,int]]:
    out=[]; out_q=[]; counts={'sub':0,'ins':0,'del':0}
    seq=clean_dna(seq)
    for i,b in enumerate(seq):
        q=int(qscores[i]) if i<len(qscores) else Q_MAX
        if rng.random()<float(deletion_rate):
            counts['del']+=1
        else:
            p_sub=_q_to_error_probability(q) if substitution_model=='Phred Q-score driven' else float(configured_substitution_rate)
            ob=b
            if rng.random()<p_sub:
                ob=rng.choice([x for x in BASES if x!=b]); counts['sub']+=1
            out.append(ob); out_q.append(q)
        # Exactly as app_ECC: insertion is independent of substitution/Q model.
        if rng.random()<float(insertion_rate):
            out.append(rng.choice(BASES)); out_q.append(q); counts['ins']+=1
    return ''.join(out),out_q,counts


def simulate_fastq(
    strands:List[Dict[str,Any]],
    coverage:int=10,
    sub_rate:float=0.003,
    ins_rate:float=0.0002,
    del_rate:float=0.0002,
    mean_q:int=30,
    seed:int=7,
    error_probability_model:str='Independent Bernoulli per nucleotide',
    qscore_profile:str='Fixed Q',
    q_mean:float=30.0,
    q_std:float=4.0,
    q_start:float=35.0,
    q_end:float=20.0,
    q_fixed:float=30.0,
    strand_dropout_rate:float=0.0,
)->Tuple[bytes,Dict]:
    """Simulate FASTQ from prepared strands without encoding a hidden SID in headers.

    For normal prepared rows, only the Payload/BODY region receives S/I/D.
    FBR, SI, filler and RBR are kept structural, matching app_ECC.  HEDges rows
    therefore remain full 125-nt reads while the 76-nt BODY receives errors.
    If only ``Full strand`` is supplied, the whole sequence is treated as payload.

    Error probability model:
      * Independent Bernoulli per nucleotide: configured substitution rate.
      * Phred Q-score driven: p_sub = 10^(-Q/10) using one of the three Q profiles.
    Insertions and deletions remain independent Bernoulli events in both modes.
    """
    coverage=max(1,int(coverage))
    model=str(error_probability_model or 'Independent Bernoulli per nucleotide')
    if model not in {'Independent Bernoulli per nucleotide','Phred Q-score driven'}:
        raise ValueError(f'Unsupported error probability model: {model}')
    rng=random.Random(int(seed))
    lines=[]; total={'records':0,'sub':0,'ins':0,'del':0,'dropped':0}; per=[]
    global_read=0; all_q=[]
    for fallback,row in enumerate(strands):
        sid=int(row.get('Strand ID', fallback))
        fbr=clean_dna(row.get('FBR','')); si=clean_dna(row.get('SI',''))
        payload=clean_dna(row.get('Payload','')); filler=clean_dna(row.get('Filler','')); rbr=clean_dna(row.get('RBR',''))
        has_segments=bool(fbr or si or rbr or ('Payload' in row))
        if not has_segments:
            payload=clean_dna(row.get('Full strand',''))
            fbr=si=filler=rbr=''
        if rng.random()<float(strand_dropout_rate):
            total['dropped']+=1
            per.append({'Strand ID':sid,'Reads':0,'Substitutions':0,'Insertions':0,'Deletions':0,'Dropped':'Yes'})
            continue
        sc={'Strand ID':sid,'Reads':0,'Substitutions':0,'Insertions':0,'Deletions':0,'Dropped':'No'}
        for copy_no in range(1,coverage+1):
            read_rng=random.Random(f'{seed}|strand={sid}|copy={copy_no}')
            if model=='Phred Q-score driven':
                pq=make_qscores(len(payload),read_rng,qscore_profile,q_mean,q_std,q_start,q_end,q_fixed)
            else:
                pq=[_p_to_q(float(sub_rate))]*len(payload)
            noisy_payload,noisy_q,c=_mutate_segment(payload,read_rng,model,float(sub_rate),pq,float(ins_rate),float(del_rate))
            full_read=fbr+si+noisy_payload+filler+rbr
            full_q=[STRUCTURAL_Q]*(len(fbr)+len(si))+noisy_q+[STRUCTURAL_Q]*(len(filler)+len(rbr))
            if len(full_q)!=len(full_read): raise RuntimeError('FASTQ quality length mismatch')
            global_read+=1
            # Generic read header deliberately contains no strand ID. Decoder must use sequence SI,
            # or an ID that is explicitly present in an externally supplied FASTQ header.
            lines += [f'@read{global_read}|copy{copy_no}', full_read, '+', ''.join(_qchar(q) for q in full_q)]
            total['records']+=1
            for k in ('sub','ins','del'): total[k]+=c[k]
            sc['Reads']+=1; sc['Substitutions']+=c['sub']; sc['Insertions']+=c['ins']; sc['Deletions']+=c['del']
            all_q.extend(pq)
        per.append(sc)
    stats={
        'summary':total,'per_strand':per,'coverage':coverage,
        'error_probability_model':model,
        'qscore_profile':qscore_profile if model=='Phred Q-score driven' else 'Fixed from Bernoulli p',
        'configured_substitution_rate':float(sub_rate),
        'configured_insertion_rate':float(ins_rate),
        'configured_deletion_rate':float(del_rate),
        'configured_dropout_rate':float(strand_dropout_rate),
        'mean_payload_q':(sum(all_q)/len(all_q)) if all_q else 0.0,
        'min_payload_q':min(all_q) if all_q else 0,
        'max_payload_q':max(all_q) if all_q else 0,
        'q_mean':float(q_mean),'q_std':float(q_std),'q_start':float(q_start),'q_end':float(q_end),'q_fixed':float(q_fixed),
        'assignment_signal':'Embedded SI in sequence; FASTQ headers contain no simulator strand ID',
    }
    return ('\n'.join(lines)+'\n').encode('utf-8'), stats

def _nw_align(a:str,b:str)->Tuple[str,str,List[int|None]]:
    """Global alignment; returns aligned a, aligned b and b-index per aligned column."""
    n,m=len(a),len(b); match=2; mismatch=-1; gap=-2
    dp=[[0]*(m+1) for _ in range(n+1)]
    bt=[[None]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=i*gap; bt[i][0]='U'
    for j in range(1,m+1): dp[0][j]=j*gap; bt[0][j]='L'
    for i in range(1,n+1):
        ai=a[i-1]
        for j in range(1,m+1):
            diag=dp[i-1][j-1]+(match if ai==b[j-1] else mismatch)
            up=dp[i-1][j]+gap
            left=dp[i][j-1]+gap
            best=max(diag,up,left)
            dp[i][j]=best
            bt[i][j]='D' if best==diag else ('U' if best==up else 'L')
    aa=[];bb=[];bi=[];i=n;j=m
    while i or j:
        t=bt[i][j]
        if i and j and t=='D': aa.append(a[i-1]); bb.append(b[j-1]); bi.append(j-1); i-=1;j-=1
        elif i and (j==0 or t=='U'): aa.append(a[i-1]); bb.append('-'); bi.append(None); i-=1
        else: aa.append('-'); bb.append(b[j-1]); bi.append(j-1); j-=1
    return ''.join(reversed(aa)),''.join(reversed(bb)),list(reversed(bi))


def _project_to_anchor(anchor:str, seq:str, q:List[int]):
    aa,bb,bi=_nw_align(anchor,seq)
    # anchor positions and insertion slots 0..len(anchor)
    bases=[None]*len(anchor); bq=[None]*len(anchor); inserts=defaultdict(list)
    apos=0
    for ca,cb,idx in zip(aa,bb,bi):
        if ca=='-':
            slot=apos
            if cb!='-': inserts[slot].append((cb, q[idx] if idx is not None and idx<len(q) else 0))
        else:
            if cb!='-':
                bases[apos]=cb; bq[apos]=q[idx] if idx is not None and idx<len(q) else 0
            apos+=1
    return bases,bq,inserts


def _phred_error(q:float)->float:
    return 10.0**(-float(q)/10.0)


def _vote_base(obs:List[Tuple[str,int]], method:str)->Tuple[str,float,bool,Dict[str,float]]:
    if not obs: return 'N',0.0,False,{b:0.0 for b in BASES}
    counts=Counter(b for b,_ in obs if b in BASES)
    empirical={b:counts.get(b,0)/len(obs) for b in BASES}
    if method=='Q-score likelihood voting':
        logs={}
        for cand in BASES:
            s=0.0
            for b,q in obs:
                pe=min(0.999999,max(1e-12,_phred_error(q)))
                p=(1-pe) if b==cand else pe/3.0
                s+=math.log(max(p,1e-300))
            logs[cand]=s
        mx=max(logs.values()); ex={b:math.exp(logs[b]-mx) for b in BASES}; z=sum(ex.values()) or 1
        score={b:ex[b]/z for b in BASES}
    else:
        score=empirical
    best=max(score.values()); winners=[b for b in BASES if abs(score[b]-best)<1e-12]
    return winners[0],best,len(winners)>1,score


def build_consensus(records:List[Dict[str,Any]], method:str='Q-score likelihood voting')->Tuple[str,Dict,List[Dict]]:
    if not records: raise ValueError('No reads for consensus')
    lengths=[len(r['sequence']) for r in records]
    modal=Counter(lengths).most_common(1)[0][0]
    anchor_rec=max(records,key=lambda r:(-abs(len(r['sequence'])-modal), sum(r['phred'])/max(1,len(r['phred']))))
    anchor=anchor_rec['sequence']
    projected=[_project_to_anchor(anchor,r['sequence'],r['phred']) for r in records]
    out=[]; posrows=[]; ties=0; low=0
    # insertion before each anchor base, and after final base
    for slot in range(len(anchor)+1):
        # Collect first insertion at slot from each read. Multiple-base insertions are uncommon; support columns iteratively.
        slot_lists=[]
        for bases,bq,ins in projected:
            slot_lists.append(ins.get(slot,[]))
        maxins=max([len(x) for x in slot_lists]+[0])
        for k in range(maxins):
            obs=[lst[k] for lst in slot_lists if k<len(lst)]
            # require insertion support from >50% reads to avoid preserving sporadic sequencing insertions
            if len(obs) > len(records)/2:
                b,p,t,_=_vote_base(obs,method); out.append(b); ties+=int(t); low+=int(p<0.70)
        if slot==len(anchor): break
        obs=[]; gaps=0
        for bases,bq,ins in projected:
            if bases[slot] is None: gaps+=1
            else: obs.append((bases[slot], int(bq[slot] or 0)))
        if gaps > len(records)/2:
            continue
        b,p,t,scores=_vote_base(obs,method)
        if b!='N': out.append(b)
        ties+=int(t); low+=int(p<0.70)
        posrows.append({'Position':slot+1,'Coverage':len(obs),'Gap votes':gaps,'Consensus':b,'Confidence':p,'Tie':t,
                        'A':sum(x[0]=='A' for x in obs),'C':sum(x[0]=='C' for x in obs),'G':sum(x[0]=='G' for x in obs),'T':sum(x[0]=='T' for x in obs)})
    seq=''.join(out)
    return seq, {'reads':len(records),'anchor_length':len(anchor),'consensus_length':len(seq),'tie_positions':ties,'low_confidence_positions':low,'method':method}, posrows


def consensus_from_fastq(data:bytes|str, expected_sids:List[int]|None=None, assignment_mode:str='Header strand ID', reads_per_strand:int=10, method:str='Q-score likelihood voting'):
    # Compatibility helper for external FASTQ that explicitly stores strand IDs in headers.
    # Record-order assignment is intentionally disabled in V4 because it is a hidden positional hint.
    if assignment_mode != 'Header strand ID':
        raise ValueError('Record-order FASTQ assignment is disabled; use embedded SI or an explicit FASTQ header strand ID')
    recs=parse_fastq(data)
    groups=defaultdict(list); unassigned=0
    for i,r in enumerate(recs):
        sid=sid_from_header(r['header'])
        if sid is None: unassigned+=1; continue
        groups[int(sid)].append(r)
    results=[]; consensus={}; positions=[]
    for sid in sorted(groups):
        seq,stats,pos=build_consensus(groups[sid],method)
        consensus[sid]=seq
        results.append({'Strand ID':sid,'Coverage':stats['reads'],'Consensus length':len(seq),'Ties':stats['tie_positions'],'Low-confidence':stats['low_confidence_positions'],'Method':method})
        for row in pos: positions.append({'Strand ID':sid,**row})
    return consensus,results,positions,{'records':len(recs),'usable_reads':sum(len(v) for v in groups.values()),'unassigned_reads':unassigned,'covered_strands':len(groups)}


def consensus_from_intrinsic_ids(
    data: bytes | str,
    *,
    fbr: str,
    rbr: str,
    si_nt: int,
    si_decoder,
    method: str = 'Q-score likelihood voting',
    target: str = 'body',
):
    """Build FASTQ consensus without record-order or session/reference assignment.

    Read assignment priority:
      1) SI physically present in a full prepared read (intrinsic sequence data),
      2) explicit strand ID physically present in an external FASTQ header.

    ``target='body'`` strips FBR|SI and RBR before voting (HEDges).
    ``target='full'`` votes on the full read (legacy SM/rule-base strand layer).
    """
    if target not in {'body','full'}:
        raise ValueError("target must be 'body' or 'full'")
    recs = parse_fastq(data)
    fbr = clean_dna(fbr); rbr = clean_dna(rbr); si_nt = int(si_nt)
    groups = defaultdict(list); sources = Counter(); unassigned = 0
    for rec in recs:
        seq = clean_dna(rec.get('sequence','')); q = list(rec.get('phred') or [])
        sid = None; target_rec = rec; source = None
        min_full = len(fbr) + si_nt + len(rbr)
        if seq.startswith(fbr) and seq.endswith(rbr) and len(seq) >= min_full:
            obs_si = seq[len(fbr):len(fbr)+si_nt]
            try:
                sid = int(si_decoder(obs_si))
            except Exception:
                sid = None
            if sid is not None:
                source = 'Embedded SI'
                if target == 'body':
                    a = len(fbr) + si_nt
                    b = len(seq) - len(rbr) if rbr else len(seq)
                    target_rec = {**rec, 'sequence': seq[a:b], 'phred': q[a:b]}
        if sid is None:
            header_sid = sid_from_header(rec.get('header',''))
            if header_sid is not None:
                sid = int(header_sid); source = 'FASTQ header'
                target_rec = rec
        if sid is None:
            unassigned += 1
            continue
        groups[sid].append(target_rec); sources[source] += 1

    consensus = {}; results = []; positions = []
    for sid in sorted(groups):
        seq, stats, pos = build_consensus(groups[sid], method)
        consensus[sid] = seq
        results.append({
            'Strand ID': sid, 'Coverage': stats['reads'], 'Consensus length': len(seq),
            'Ties': stats['tie_positions'], 'Low-confidence': stats['low_confidence_positions'],
            'Method': method,
        })
        for row in pos:
            positions.append({'Strand ID': sid, **row})
    return consensus, results, positions, {
        'records': len(recs),
        'usable_reads': sum(len(v) for v in groups.values()),
        'unassigned_reads': unassigned,
        'covered_strands': len(groups),
        'assignment': 'Embedded SI → explicit FASTQ header fallback',
        'assignment_counts': dict(sources),
    }
