"""Reproduce one illustrative run, not a claim of population-level accuracy."""
from dataclasses import asdict
from pathlib import Path
import csv
import json
import hedges_bridge as hb
import process_model as pm
import read_pipeline as rp
import evaluation as ev


def main():
    out=Path('visual_demo_output');out.mkdir(exist_ok=False)
    data=b'Visual DNA storage demo: inspect reads, consensus, and HEDGES.'
    strands,_=hb.encode_bytes(data);rows=[r.row() for r in strands]
    cfg=pm.Settings(seed=17,founders=5,coverage=5,seq_sub=.02,seq_ins=.002,seq_del=.002)
    sim=pm.simulate(rows,cfg)
    public=rp.Settings(stored_bytes=len(data),expected_strands=len(rows),max_edits=1)
    result=rp.run(sim['r1'],public,decode_file=False)
    reads=ev.sequencing(rows,sim);consensus=ev.consensus_report(rows,result)
    result.update(rp.reconstruct(result['consensus'],public))
    summary={'purpose':'Single illustrative run; synthetic rates, no calibrated vendor claim.',
             'simulation_settings':asdict(cfg),'strands':len(rows),
             'read_edits':reads['totals'],'consensus':consensus,
             'hedges':ev.hedges_report(rows,result),'file_exact':result['data']==data}
    (out/'R1.fastq').write_bytes(sim['r1'])
    (out/'decoder_settings.json').write_text(json.dumps(asdict(public),indent=2))
    (out/'measurements.json').write_text(json.dumps(summary,indent=2))
    with (out/'read_errors.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=reads['reads'][0]);w.writeheader();w.writerows(reads['reads'])
    if result['data'] is not None:(out/'recovered.txt').write_bytes(result['data'])
    print(json.dumps({'strands':len(rows),'read_edits':reads['totals'],
                      'before_rate':consensus['before_rate'],'after_rate':consensus['after_rate'],
                      'reduction_percent':consensus['reduction_percent'],'file_exact':summary['file_exact']},indent=2))


if __name__=='__main__':main()
