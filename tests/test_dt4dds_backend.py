from dataclasses import replace

import pytest

import dt4dds_backend as db
import hedges_bridge as hb


pytestmark=pytest.mark.skipif(not db.available(),reason='DT4DDS 1.1.0 required')


@pytest.fixture
def rows():
    return [s.row() for s in hb.encode_bytes(b'DNA test')[0]]


def test_official_backend_reproducible_and_abundance_aware(rows):
    cfg=replace(db.MolecularSettings(),synthesis_coverage=20,pcr_cycles=3,aliquot_molecules_per_strand=10)
    a=db.run(rows,cfg);b=db.run(rows,cfg)
    assert a['_templates']==b['_templates']
    assert a['stages']==b['stages']
    assert [x['Stage'] for x in a['stages']]==['Synthesis','Storage','PCR','Sampling']
    assert a['stages'][0]['molecules']>0
    assert a['stages'][1]['molecules']==a['stages'][0]['molecules']
    assert a['stages'][2]['molecules']>a['stages'][1]['molecules']
    assert a['stages'][3]['molecules']==10
    assert sum(x['copies'] for x in a['_templates'])==10
    assert all(x['measurement'].startswith('abundance-weighted') for x in a['stages'])


def test_manifest_pins_sources_and_discloses_measurement(rows):
    run=db.run(rows,replace(db.MolecularSettings(),synthesis_coverage=10,pcr_cycles=2,aliquot_molecules_per_strand=5))
    manifest=run['preset']
    assert manifest['molecular_engine']=='dt4dds'
    assert manifest['dt4dds_version']=='1.1.0'
    assert set(manifest['source_sha256'])=={'array_synthesis_twist.yaml','aging.yaml','pcr_Q5.yaml','properties.yaml'}
    assert all(len(v)==64 for v in manifest['source_sha256'].values())
    assert 'not counts of newly generated events' in manifest['measurement_note']


def test_no_simplified_fallback(monkeypatch):
    def fail():
        raise RuntimeError('DT4DDS unavailable')
    monkeypatch.setattr(db,'_imports',fail)
    assert db.available() is False
