"""Post-hoc measurement only. Never imported by read_pipeline or the decoder.

Edit counts are a deterministic minimum-edit alignment, not event provenance.
Ties prefer diagonal, deletion, insertion. N is unresolved, counted separately.
"""
from collections import Counter
from functools import lru_cache
import read_pipeline as rp
import process_model as pm
import hedges_bridge as hb
from rapidfuzz.distance import Levenshtein


@lru_cache(maxsize=30000)
def alignment(reference, observed, partial=False):
    n,m=len(reference),len(observed)
    d=[[0]*(m+1) for _ in range(n+1)]
    for i in range(n+1):d[i][0]=0 if partial else i
    for j in range(m+1):d[0][j]=j
    for i in range(1,n+1):
        for j in range(1,m+1):
            d[i][j]=min(d[i-1][j-1]+(reference[i-1]!=observed[j-1]),d[i-1][j]+1,d[i][j-1]+1)
    i=min(range(n+1),key=lambda k:(d[k][m],-k)) if partial else n
    j=m;a=[];b=[]
    while j or (i and not partial):
        if i and j and d[i][j]==d[i-1][j-1]+(reference[i-1]!=observed[j-1]):
            a.append(reference[i-1]);b.append(observed[j-1]);i-=1;j-=1
        elif i and d[i][j]==d[i-1][j]+1:
            a.append(reference[i-1]);b.append('-');i-=1
        else:a.append('-');b.append(observed[j-1]);j-=1
    return ''.join(a[::-1]),''.join(b[::-1])


@lru_cache(maxsize=16384)
def edits(reference, observed, partial=False):
    if reference==observed:
        return {'S':0,'I':0,'D':0,'N':observed.count('N'),'reference_nt':len(reference),'errors':observed.count('N'),'rate':observed.count('N')/len(reference) if reference else None}
    # Most ART reads use a full-reference comparison.  RapidFuzz performs this
    # global Levenshtein evaluation in compiled code, which keeps the complete
    # pre-consensus measurement practical at normal coverage.  The explicit
    # dynamic-programming alignment remains necessary for semi-global reads.
    if not partial:
        c={'S':0,'I':0,'D':0,'N':0,'reference_nt':len(reference)}
        for op in Levenshtein.editops(reference,observed):
            if op.tag=='replace':
                c['N' if observed[op.dest_pos]=='N' else 'S']+=1
            elif op.tag=='insert':
                c['N' if observed[op.dest_pos]=='N' else 'I']+=1
            elif op.tag=='delete':c['D']+=1
        c['errors']=sum(c[k] for k in ('S','I','D','N'))
        c['rate']=c['errors']/c['reference_nt'] if c['reference_nt'] else None
        return c
    a,b=alignment(reference,observed,partial)
    c={'S':0,'I':0,'D':0,'N':0,'reference_nt':len(a.replace('-',''))}
    for x,y in zip(a,b):
        if y=='N':c['N']+=1
        elif x=='-':c['I']+=1
        elif y=='-':c['D']+=1
        elif x!=y:c['S']+=1
    c['errors']=sum(c[k] for k in ('S','I','D','N'))
    c['rate']=c['errors']/c['reference_nt'] if c['reference_nt'] else None
    return c


def strict_similarity(reference,observed):
    """User-defined positional score: unequal lengths receive zero; no alignment."""
    if not reference or len(reference)!=len(observed):return 0.0
    return 100*sum(a==b for a,b in zip(reference,observed))/len(reference)


def sequencing(rows,run):
    refs={int(r['Strand ID']):r['Full strand'] for r in rows}
    mates=[rp.parse(run['r1'])]+([rp.parse(run['r2'])] if run['r2'] else [])
    records=[]
    for truth in run['truth']:
        for mate,reads in enumerate(mates):
            rec=reads[truth['unit']];ref=refs[truth['sid']]
            if mate:ref=pm.rc(ref)
            partial=run['settings']['read_length']<len(truth['template'])
            measures=edits(ref,rec['seq'],partial)
            molecular=truth.get('sequencing_edits',[])
            sequencing_counts=molecular[mate] if molecular else None
            records.append({'unit':truth['unit'],'Strand ID':truth['sid'],'mate':mate+1,
                            'read_name':rec['name'],'length':len(rec['seq']),
                            'mean_Q':sum(rec['q'])/len(rec['q']),
                            'Strict similarity (%)':strict_similarity(ref[:run['settings']['read_length']],rec['seq']),
                            'Expected read nt':len(ref[:run['settings']['read_length']]),
                            'comparison':'Observed span (semi-global)' if partial else 'Full strand',
                            **measures,**({'Sequencing S':sequencing_counts['S'],'Sequencing I':sequencing_counts['I'],'Sequencing D':sequencing_counts['D']} if sequencing_counts else {}),
                            'sequence':rec['seq'],'reference':ref})
    units=Counter(t['sid'] for t in run['truth'])
    read_bases=Counter()
    for r in records:read_bases[r['Strand ID']]+=r['length']
    coverage=[{'Strand ID':sid,'Read units':units[sid],
               'FASTQ records':units[sid]*len(mates),
               'Base coverage':read_bases[sid]/len(ref)} for sid,ref in refs.items()]
    return {'reads':records,'coverage':coverage,
            'totals':{k:sum(r[k] for r in records) for k in ('S','I','D','N','errors','reference_nt')}}


def consensus_report(rows,result):
    refs={int(r['Strand ID']):str(r['Payload'])+str(r.get('Filler','')) for r in rows}
    all_bodies=result.get('all_consensus',result['consensus']);table=[]
    for sid,ref in refs.items():
        reads=[r['sequence'] for r in result.get('observed_bodies',[]) if r['sid']==sid]
        seq=all_bodies.get(sid,all_bodies.get(str(sid),''))
        before=[edits(ref,s) for s in reads]
        after=edits(ref,seq) if seq else None
        raw=sum(e['errors'] for e in before);den=len(ref)*len(reads)
        row={'Strand ID':sid,'Reads assigned':len(reads),'Consensus':seq,
             'Outcome':'NO CONSENSUS' if not seq else 'UNCERTAIN' if 'N' in seq else 'EXACT' if seq==ref else 'ERRORS REMAIN',
             'Before errors (sum)':raw,'Before reference nt':den,
             'Before rate':raw/den if den else None,
             'After rate':after['rate'] if after else None}
        for k in ('S','I','D','N'):row['Before '+k]=sum(e[k] for e in before);row['After '+k]=after[k] if after else None
        row['Reduction %']=100*(1-row['After rate']/row['Before rate']) if row['Before rate'] and row['After rate'] is not None else None
        table.append(row)
    matched=[r for r in table if r['Before rate'] is not None and r['After rate'] is not None]
    before=sum(r['Before rate'] for r in matched)/len(matched) if matched else None
    after=sum(r['After rate'] for r in matched)/len(matched) if matched else None
    return {'strands':table,'matched_strands':len(matched),'total_strands':len(refs),
            'before_rate':before,'after_rate':after,
            'reduction_percent':100*(1-after/before) if before and after is not None else None,
            'definition':'Mean per-strand edit rate over the same nonempty consensus groups. Before: all assigned bodies vs the assigned strand reference, including outliers/capped reads. After: one consensus per group. N counts as unresolved error. Missing groups are excluded from rates and reported separately.'}


@lru_cache(maxsize=2048)
def clean_payload(body,sid,mode):
    decoded=hb.decode_body(body,sid,0,mode)
    return decoded.payload_bits if decoded.status=='OK' else None


def hedges_report(rows,result):
    if result['settings']['codec']!='TVH-V5.2 HEDges':return []
    calls={int(r['strand_id']):r for r in result.get('strand_results',[])}
    table=[]
    for row in rows:
        sid=int(row['Strand ID']);call=calls.get(sid,{})
        bits=call.get('payload_bits');truth=clean_payload(row['Payload'],sid,result['settings']['sid_mode'])
        exact=bits==truth if bits is not None and truth is not None else None
        table.append({'Strand ID':sid,'Decoder status':call.get('status','NOT DECODED'),
                      'Payload comparison':'EXACT' if exact else 'WRONG' if exact is False else 'UNAVAILABLE',
                      'Wrong payload bits':sum(x!=y for x,y in zip(bits,truth))+abs(len(bits)-len(truth)) if exact is not None else None,
                      'Decoder edit distance':call.get('edit_distance')})
    return table
