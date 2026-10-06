"""Produce a small example with REAL ART; fails clearly when ART is absent."""
from pathlib import Path
from dataclasses import asdict
import io,json,zipfile,hashlib
import pandas as pd
import hedges_bridge as hb
import read_pipeline as rp
import art_sequencing as art
import channel_runs as cr


def main():
    data=b'DNA test'
    strands,_=hb.encode_bytes(data);rows=[s.row() for s in strands]
    settings=rp.Settings(expected_strands=len(rows),stored_bytes=len(data))
    run=cr.finish(rows,art.ArtSettings(),reads_per_oligo=20)
    result=rp.run(run['r1'],settings,run['r2'])
    if result['data']!=data:raise RuntimeError('Example did not recover exactly')
    dest=Path(__file__).parent/'art_example';dest.mkdir(exist_ok=True)
    files={'R1.fastq':run['r1'],'R2.fastq':run['r2'],
        'decoder_settings.json':json.dumps(asdict(settings),indent=2).encode(),
        'evaluation/reference_strands.csv':pd.DataFrame(rows).to_csv(index=False).encode(),
        'evaluation/expected_hash.txt':hashlib.sha256(data).hexdigest().encode(),
        'evaluation/channel_config.json':json.dumps({'art':{k:v for k,v in run['art'].items() if k!='templates'}},indent=2).encode()}
    with zipfile.ZipFile(dest/'sequencing_output.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name,content in files.items():z.writestr(name,content)
    (dest/'prepared_strands.csv').write_bytes(files['evaluation/reference_strands.csv'])
    (dest/'decoder_settings.json').write_bytes(files['decoder_settings.json'])
    (dest/'expected_stored_sha256.txt').write_bytes(files['evaluation/expected_hash.txt'])
    (dest/'README.txt').write_text('Real ART HS25, PE125, 20 pairs/oligo, master seed 12345. Input bytes: DNA test.\nUpload sequencing_output.zip on Page 2; apply public settings; Build consensus; Decode.\nExpected stored-byte SHA-256: '+hashlib.sha256(data).hexdigest()+'\nThis smoke test is not a reliability estimate.\n')
    print('Real ART example: exact stored-byte recovery; '+str(dest/'sequencing_output.zip'))

if __name__=='__main__':main()
