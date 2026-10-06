import os
from dataclasses import replace
from collections import Counter
import pytest
import art_sequencing as art
import channel_runs as cr
import channel_presets as cp
import process_model as pm
import read_pipeline as rp
import hedges_bridge as hb

@pytest.fixture
def rows():
    strands,_=hb.encode_bytes(b'DNA test')
    return [s.row() for s in strands]

@pytest.fixture
def real_art():
    if not art.detect_art_executable():pytest.skip('Real ART executable not available')
    return art.ArtSettings()

@pytest.mark.parametrize('paired',[False,True])
def test_real_art_counts_reproducibility_and_decode(rows,real_art,paired):
    cfg=replace(real_art,paired=paired)
    a=cr.finish(rows,cfg,reads_per_oligo=20);b=cr.finish(rows,cfg,reads_per_oligo=20)
    assert a['r1']==b['r1'] and a['r2']==b['r2']
    assert len(rp.parse(a['r1']))==20
    assert len(rp.parse(a['r2']))==(20 if paired else 0)
    assert all(x['name'].startswith('read') and 'v000' not in x['name'] for x in rp.parse(a['r1']))
    assert a['art']['logs'] and 'ART' in a['art']['logs'][0]['stdout']
    out=rp.run(a['r1'],rp.Settings(expected_strands=1,stored_bytes=8),a['r2'])
    assert out['data']==b'DNA test'
    assert len(a['truth'])==20


def test_art_missing_does_not_fabricate(rows):
    with pytest.raises(FileNotFoundError,match='not found'):
        art.run_art_illumina(rows,art_exe='/definitely_missing/art_illumina')


def test_sampling_preserves_variants_and_bias():
    pool=[{'sid':0,'template':'A'*125,'copies':900}, {'sid':0,'template':'C'*125,'copies':90}, {'sid':1,'template':'G'*125,'copies':10}]
    out=art.allocate(pool,10000,71)
    assert sum(x['count'] for x in out)==10000
    assert out[0]['count']>8500 and out[2]['count']<200
    assert out==art.allocate(pool,10000,71)
    equal=art.allocate(pool,10000,71,'Equal per strand')
    by_sid=Counter()
    for row in equal:by_sid[row['sid']]+=row['count']
    assert by_sid=={0:5000,1:5000}
    major=next(row for row in equal if row['sid']==0 and row['template']=='A'*125)['count']
    minor=next(row for row in equal if row['sid']==0 and row['template']=='C'*125)['count']
    assert major>minor


def test_real_short_templates_are_never_dropped(real_art):
    pool=[{'sid':0,'template':'ACGT'*30+'ACG','count':7}]
    with pytest.raises(ValueError,match='No templates were removed'):
        art.run_art(pool,replace(real_art,short_policy='Stop if template is short'))
    out=art.run_art(pool,real_art)
    assert len(rp.parse(out['r1']))==7
    assert {len(x['seq']) for x in rp.parse(out['r1'])}=={123}
    assert out['art']['shortened_read_units']==7


def test_cascade_real_art_inheritance(rows,real_art):
    cfg=replace(cp.defaults(),cycles=2,syn_sub=1.,syn_del=0.,pcr_sub=0.,pcr_ins=0.,pcr_del=0.,scope='BODY only')
    checkpoint=None
    for stage in ['Synthesis','Storage','PCR']:
        checkpoint=pm.simulate(rows,cfg,stop_after=stage,checkpoint=checkpoint)
    population=art.pool_templates(rows,checkpoint)
    out=cr.finish(rows,real_art,checkpoint=checkpoint,total_units=40)
    assert len(out['truth'])==40
    assert all(t['template'] in [x['template'] for x in population] for t in out['truth'])
    assert all(t['template']!=rows[0]['Full strand'] for t in out['truth'])
    assert sum(x['Read units'] for x in out['abundance'])==40


def test_abundance_and_survival_checkpoints(rows):
    rows=[dict(rows[0],**{'Strand ID':i}) for i in range(10)]
    cfg=replace(cp.defaults(),founders=20,storage=True,storage_survival=0.,cycles=1)
    syn=pm.simulate(rows,cfg,stop_after='Synthesis')
    assert sum(syn['_initial_counts'])==200 and len(set(syn['_initial_counts']))>1
    store=pm.simulate(rows,cfg,stop_after='Storage',checkpoint=syn)
    assert all(not x for x in store['_populations'])
    pcr=pm.simulate(rows,cfg,stop_after='PCR',checkpoint=store)
    with pytest.raises(ValueError,match='No surviving'):
        art.allocate(art.pool_templates(rows,pcr),100,1)


def test_preset_manifest():
    cfg=cp.defaults();m=cp.manifest(cfg)
    assert not m['overrides'] and m['source_sha256']
    assert cfg.pcr_sub==1.09e-4/280
    assert 'Custom' in cp.manifest(replace(cfg,syn_del=.1))['label']


def test_aging_units_and_sequential_equivalence(rows):
    cfg=replace(cp.defaults(),storage=True,storage_halflives=2.,storage_rate_per_halflife=.1,cycles=1)
    whole=pm.simulate(rows,cfg,stop_after='PCR')
    step=None
    for stage in ['Synthesis','Storage','PCR']:
        step=pm.simulate(rows,cfg,stop_after=stage,checkpoint=step)
    assert whole['_populations']==step['_populations']
    assert next(r for r in step['stages'] if r['Stage']=='Storage')['sub']>0


def test_art_evaluation_counts_and_package(rows,real_art):
    import io,json,zipfile
    import imports,evaluation
    sim=cr.finish(rows,replace(real_art,insertion_r1=.002,deletion_r1=.002),reads_per_oligo=30)
    report=evaluation.sequencing(rows,sim)
    stage=sim['sequencing_stage']
    for letter,key in [('S','sub'),('I','ins'),('D','del')]:
        assert sum(r['Sequencing '+letter] for r in report['reads'])==stage[key]
    package=io.BytesIO()
    with zipfile.ZipFile(package,'w') as z:
        z.writestr('R1.fastq',sim['r1']);z.writestr('R2.fastq',sim['r2'])
        z.writestr('evaluation/channel_config.json',json.dumps({'engine':'ART'}))
    loaded=imports.package(package.getvalue())
    assert loaded['R1.fastq']==sim['r1']


def test_wilson_does_not_claim_99_from_few_successes():
    from run_channel_experiment import wilson
    lo,hi=wilson(30,30)
    assert lo<.99 and hi<=1
