"""Bounded, reference-free search of observed consensus alternatives.

Unique means unique within the declared evidence alternatives and edit radius.
A hard per-strand wall-clock budget includes spawn/startup. Incomplete search
never yields an accepted payload. Kept-base errors remain covered by HEDGES K.
"""
import itertools
import math
import multiprocessing as mp
import time
from dataclasses import asdict
import hedges_bridge as hb


def options(sequence,positions):
    if not positions:
        return [list('ACGT') if b=='N' else [b] for b in sequence]
    rebuilt=''.join(p['base'] for p in positions if p['base']!='-')
    if rebuilt!=sequence:raise ValueError('Consensus evidence does not match sequence')
    return [p.get('alternatives',list('ACGT')) if p['base']=='N' else ['' if p['base']=='-' else p['base']] for p in positions]


def search(sequence,positions,sid,c):
    opts=options(sequence,positions);uncertain=sum(len(x)>1 for x in opts)
    total=math.prod(len(x) for x in opts)
    row={'strand_id':sid,'hypotheses':total,'tested':0,'uncertain_positions':sequence.count('N'),'payload_bits':None}
    if 'N' in sequence and not c.rescue_uncertain:return {**row,'status':'UNCERTAIN'}
    if sequence.count('N')>c.max_uncertain or total>c.max_candidates:return {**row,'status':'SEARCH_LIMIT'}
    payloads={};states=0;seen=set()
    for pieces in itertools.product(*opts):
        body=''.join(pieces)
        if body in seen:continue
        seen.add(body);row['tested']+=1
        decoded=hb.decode_body(body,sid,c.max_edits,c.sid_mode);states+=decoded.states
        if decoded.status=='AMBIGUOUS':return {**row,'status':'AMBIGUOUS','states':states}
        if decoded.status=='OK':payloads[decoded.payload_bits]=decoded
        if len(payloads)>1:return {**row,'status':'AMBIGUOUS','states':states}
    if not payloads:return {**row,'status':'FAIL','states':states}
    decoded=next(iter(payloads.values()))
    return {**row,**asdict(decoded),'states':states,'search_complete':True}


def _worker(conn,sequence,positions,sid,c):
    try:conn.send(search(sequence,positions,sid,c))
    except Exception as exc:conn.send({'strand_id':sid,'status':'ERROR','error':str(exc),'payload_bits':None})
    finally:conn.close()


def bounded_search(sequence,positions,sid,c):
    # Check combinatorial budgets before starting a worker.
    opts=options(sequence,positions);total=math.prod(len(x) for x in opts)
    if ('N' in sequence and not c.rescue_uncertain) or sequence.count('N')>c.max_uncertain or total>c.max_candidates:
        return search(sequence,positions,sid,c)
    ctx=mp.get_context('spawn');parent,child=ctx.Pipe(duplex=False)
    proc=ctx.Process(target=_worker,args=(child,sequence,positions,sid,c),daemon=True)
    start=time.monotonic();proc.start();child.close()
    try:
        remaining=max(0,c.search_seconds-(time.monotonic()-start))
        if parent.poll(remaining):
            try:row=parent.recv()
            except EOFError:row={'strand_id':sid,'status':'ERROR','payload_bits':None,'error':'Worker ended unexpectedly'}
        else:row={'strand_id':sid,'status':'SEARCH_LIMIT','payload_bits':None,'reason':'Time budget exceeded','hypotheses':total}
    finally:
        if proc.is_alive():proc.terminate()
        proc.join();parent.close();proc.close()
    return {**row,'seconds':time.monotonic()-start}


def decode(result,c,progress=None,retry_only=False):
    import read_pipeline as rp
    rp.validate(c)
    if c.codec!='TVH-V5.2 HEDges':return rp.reconstruct(result['consensus'],c)
    bodies=result['all_consensus'];expected=list(range(c.expected_strands)) if c.expected_strands else list(range(max(bodies,default=-1)+1))
    verified={};rows=[]
    prior={r['strand_id']:r for r in result.get('search_results',[])}
    old=result.get('decode_settings',{})
    if old.get('max_edits')!=c.max_edits or old.get('sid_mode')!=c.sid_mode:retry_only=False
    for index,sid in enumerate(expected):
        sequence=bodies.get(sid,'')
        positions=[p for p in result['positions'] if p['Strand ID']==sid]
        if retry_only and sid in prior and prior[sid]['status']!='SEARCH_LIMIT':row=prior[sid]
        else:row=bounded_search(sequence,positions,sid,c) if sequence else {'strand_id':sid,'status':'MISSING','payload_bits':None}
        rows.append(row)
        if row['status']=='OK':
            verified[sid]=hb.HedgesDecoded(**{k:row[k] for k in hb.HedgesDecoded.__dataclass_fields__})
        if progress:progress((index+1)/len(expected))
    extras=set(bodies)-set(expected)
    if not expected or len(verified)!=len(expected) or extras:
        return {'status':'FAIL','data':None,'error':'Some strands remain unresolved. See strand results.','strand_results':rows,'search_results':rows,'decode_settings':asdict(c)}
    # Payloads are already validated by a completed search. Do not decode again
    # without a time budget or select candidates using file hash/format.
    output=rp.reconstruct({sid:bodies[sid] for sid in expected},c,verified=verified)
    output['search_results']=rows;output['strand_results']=rows;output['decode_settings']=asdict(c)
    return output
