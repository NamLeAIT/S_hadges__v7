from dataclasses import replace
import time
from types import SimpleNamespace
import numpy as np
import hedges_bridge as hb
import read_pipeline as rp
import process_model as pm
import candidate_search as cs


def test_actual_uncertain_base_and_gap_rescue():
    strands,_=hb.encode_bytes(b'DNA test');strand=strands[0];body=strand.body
    c=rp.Settings(expected_strands=1,stored_bytes=8,max_edits=0)
    sequence=body[:10]+'N'+body[11:]
    positions=[{'base':b,'alternatives':list('ACGT')} for b in sequence]
    found=cs.bounded_search(sequence,positions,0,c)
    assert found['status']=='OK' and found['payload_bits']==strand.payload_bits and found['tested']==4
    inserted=body[:10]+'N'+body[10:]
    positions=[{'base':b,'alternatives':['','A']} for b in inserted]
    found=cs.bounded_search(inserted,positions,0,c)
    assert found['status']=='OK' and found['tested']==2


def test_all_candidates_required_and_distinct_payloads_rejected(monkeypatch):
    def call(body,*args):
        return hb.HedgesDecoded(0,'OK',body,0,1,1,0)
    monkeypatch.setattr(hb,'decode_body',call)
    assert cs.search('N',[{'base':'N','alternatives':['A','G']}],0,rp.Settings())['status']=='AMBIGUOUS'
    limited=cs.search('NN',[],0,rp.Settings(max_candidates=2))
    assert limited['status']=='SEARCH_LIMIT' and limited['payload_bits'] is None


def test_hard_timeout_never_accepts_partial_search():
    # Deliberately tiny budget includes spawn/import startup, as documented.
    c=replace(rp.Settings(),search_seconds=.000001)
    start=time.monotonic();result=cs.bounded_search('A'*76,[],0,c)
    assert result['status']=='SEARCH_LIMIT' and result['payload_bits'] is None
    assert time.monotonic()-start<5


def test_pcr_bias_and_sampling_modes():
    strands,_=hb.encode_bytes(b'PCR abundance and sampling check.'*4);rows=[s.row() for s in strands]
    base=pm.Settings(seed=91,synthesis=False,sequencing=False,pcr_sub=0,cycles=8,bias_sd=.3,coverage=3)
    fixed=pm.simulate(rows,base)
    assert all(r['Read units']==3 for r in fixed['abundance'])
    assert np.std([r['PCR efficiency'] for r in fixed['abundance']])>0
    pooled=pm.simulate(rows,replace(base,sampling='Pool sampling'))
    assert sum(r['Read units'] for r in pooled['abundance'])==3*len(rows)
    assert len(set(r['Read units'] for r in pooled['abundance']))>1
    reserved=pm.simulate(rows,replace(base,sampling='Pool sampling + minimum 1'))
    assert all(r['Read units']>=1 for r in reserved['abundance'])
    assert sum(r['Read units'] for r in reserved['abundance'])==3*len(rows)
    again=pm.simulate(rows,replace(base,sampling='Pool sampling'))
    assert pooled['r1']==again['r1']
