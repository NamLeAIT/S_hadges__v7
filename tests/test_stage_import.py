from dataclasses import replace
import io,zipfile
import pytest
import hedges_bridge as hb
import process_model as pm
import read_pipeline as rp
import imports


def test_stage_error_provenance_and_quality():
    strands,_=hb.encode_bytes(b'DNA test');rows=[s.row() for s in strands]
    c=pm.Settings(synthesis=True,syn_sub=1,syn_ins=0,syn_del=0,storage=False,pcr=False,sequencing=False,founders=2,coverage=3)
    run=pm.simulate(rows,c)
    syn,storage=run['stage_quality']
    assert syn['Before errors']==0 and syn['After errors']>0
    assert syn['After exact_molecules_percent']==0
    assert storage['Before errors']==storage['After errors']==syn['After errors']
    records=rp.parse(run['r1'])
    assert all(len(r['q'])==len(r['seq']) and min(r['q'])==60 for r in records)
    clean=pm.simulate(rows,replace(c,synthesis=False,storage=True,storage_sub=1,storage_ct=0))
    assert clean['stage_quality'][0]['After errors']==0
    assert clean['stage_quality'][1]['After errors']>0


def test_package_validation():
    def blob(files):
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w') as z:
            for name,data in files.items():z.writestr(name,data)
        return b.getvalue()
    assert imports.package(blob({'R1.fastq':b'example'}))['R1.fastq']==b'example'
    with pytest.raises(ValueError):imports.package(blob({'../R1.fastq':b'bad'}))
    with pytest.raises(ValueError):imports.references(b'Strand ID,Payload\n0,ACGT\n0,AAAA\n')
