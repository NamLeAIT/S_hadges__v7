"""CLI coverage/seed sweep. Run on Windows/Linux with real ART installed.

References are used after consensus/decode only. An operational exception is
recorded separately and excluded from the recovery denominator.
"""
import argparse
from dataclasses import asdict,replace
from pathlib import Path
import json
import hashlib
import math
import time
import pandas as pd
import art_sequencing as art
import channel_runs as cr
import dt4dds_backend as db
import read_pipeline as rp
import evaluation as ev


def wilson(k,n,z=1.959963984540054):
    if not n:return None,None
    p=k/n;d=1+z*z/n
    mid=(p+z*z/(2*n))/d
    delta=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return max(0.,mid-delta),min(1.,mid+delta)


def run(args):
    from imports import references
    rows=references(Path(args.strands).read_bytes())
    if not all('Full strand' in r for r in rows):raise ValueError('Prepared strands CSV must include Full strand')
    public=rp.load_settings(Path(args.settings).read_bytes())
    if public.expected_strands!=len(rows):raise ValueError('Public expected_strands must match prepared strands for this experiment')
    if not art.detect_art_executable(args.art):raise FileNotFoundError('ART executable not found')
    expected=args.expected_hash.lower().strip()
    if len(expected)!=64 or set(expected)-set('0123456789abcdef'):raise ValueError('Provide expected stored-file SHA-256, not the pre-compression hash')
    coverage=[int(x) for x in args.coverage.split(',')]
    if any(c<1 or c>10000 or c*len(rows)>100000 for c in coverage):raise ValueError('Coverage/read budget outside supported limits')
    if args.repeats<1:raise ValueError('repeats must be positive')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    base_art=art.ArtSettings(executable=args.art,paired=args.paired,read_length=args.length)
    manifest={'arguments':vars(args),'public_settings':asdict(public),'art_settings':asdict(base_art),
              'molecular_profile':db.manifest(db.MolecularSettings(seed=args.seed)),
              'strands_sha256':hashlib.sha256(Path(args.strands).read_bytes()).hexdigest(),
              'note':'Experimental mode uses DT4DDS then total_units=coverage*designed_strands with unequal SID coverage. SE vs PE differ in nucleotide budget. All-zero-error successes do not establish a 99% reliability guarantee.'}
    (out/'experiment.json').write_text(json.dumps(manifest,indent=2))
    trials=[]
    for repetition in range(args.repeats):
        seed=args.seed+repetition
        checkpoint=None
        if args.mode=='experimental':checkpoint=db.run(rows,db.MolecularSettings(seed=seed))
        for c in coverage:
            start=time.monotonic();row={'seed':seed,'coverage':c,'paired':args.paired,'mode':args.mode}
            try:
                cfg=replace(base_art,seed=seed)
                sim=cr.finish(rows,cfg,reads_per_oligo=c,checkpoint=checkpoint,total_units=c*len(rows))
                decoded=rp.run(sim['r1'],public,sim['r2'])
                measured=ev.consensus_report(rows,decoded)
                row.update(operational_error='',file_exact=int(decoded.get('sha256')==expected),
                    consensus_exact=sum(x['Outcome']=='EXACT' for x in measured['strands']),
                    missing_consensus=sum(x['Outcome']=='NO CONSENSUS' for x in measured['strands']),
                    before_rate=measured['before_rate'],after_rate=measured['after_rate'],
                    read_units=len(sim['truth']),sequenced_nt=sum(len(x['seq']) for raw in (sim['r1'],sim['r2']) for x in rp.parse(raw)),
                    capped_units=sum(x.get('capped',0) for x in decoded['consensus_rows']),
                    decoder_status=decoded['status'])
                for status in ('OK','FAIL','AMBIGUOUS','SEARCH_LIMIT'):
                    row[status]=sum(x.get('status')==status for x in decoded.get('search_results',[]))
                artifact=out/f'art_seed{seed}_coverage{c}.json'
                artifact.write_text(json.dumps({k:v for k,v in sim['art'].items() if k!='templates'},indent=2))
            except Exception as exc:row.update(operational_error=str(exc),file_exact=None)
            row['seconds']=time.monotonic()-start;trials.append(row)
            pd.DataFrame(trials).to_csv(out/'trials.csv',index=False)
            print(row,flush=True)
    summary=[]
    for c in coverage:
        group=[r for r in trials if r['coverage']==c];good=[r for r in group if not r['operational_error']]
        n=len(good);k=sum(r['file_exact'] for r in good);lo,hi=wilson(k,n)
        summary.append({'coverage':c,'completed_trials':n,'operational_errors':len(group)-n,'exact_files':k,
                        'recovery_probability':k/n if n else None,'wilson95_low':lo,'wilson95_high':hi})
    pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--strands',required=True);p.add_argument('--settings',required=True)
    p.add_argument('--expected-hash',required=True)
    p.add_argument('--mode',choices=['experimental','illumina'],default='experimental')
    p.add_argument('--coverage',default='5,10,20,50,100,350');p.add_argument('--repeats',type=int,default=30)
    p.add_argument('--seed',type=int,default=12345);p.add_argument('--paired',action=argparse.BooleanOptionalAction,default=True)
    p.add_argument('--length',type=int,default=125);p.add_argument('--art',default='');p.add_argument('--output',default='experiment_output')
    run(p.parse_args())

if __name__=='__main__':main()
