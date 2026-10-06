"""Create a small reproducible example; ground truth remains evaluation-only."""
from pathlib import Path
from dataclasses import asdict
import json
import hashlib
import hedges_bridge as hb
import process_model as pm
import read_pipeline as rp


def main():
    target=Path('example')
    if target.exists():raise SystemExit('example already exists; preserve it or choose another working directory')
    target.mkdir()
    data=b'DNA storage example.'
    strands,meta=hb.encode_bytes(data);rows=[s.row() for s in strands]
    sim=pm.simulate(rows,pm.Settings(coverage=15,paired=True))
    (target/'R1.fastq').write_bytes(sim['r1']);(target/'R2.fastq').write_bytes(sim['r2'])
    c=rp.Settings(expected_strands=len(rows),stored_bytes=len(data))
    (target/'decoder_settings.json').write_text(json.dumps(asdict(c),indent=2))
    (target/'expected_stored_sha256.txt').write_text(hashlib.sha256(data).hexdigest()+'\n')
    (target/'simulation_evaluation.json').write_text(json.dumps({k:v for k,v in sim.items() if k not in ('r1','r2')},indent=2))
    print('Example created. Decode it with decode_fastq.py in a separate process.')


if __name__=='__main__':main()
