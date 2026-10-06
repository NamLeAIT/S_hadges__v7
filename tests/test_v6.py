from dataclasses import replace
import hedges_bridge as hb
import read_pipeline as rp
import legacy_codec as lc
import process_model as pm


def test_expanded_si_capacity_constraints_and_boundary_decode():
    mode=hb.SID_MODE_NO_HAMMING
    book=hb.codebook(hb.strand_capacity(mode),mode)
    assert len(book)==8192 and len(set(book.values()))==8192
    assert all(hb._si_bio_ok(si) for si in book.values())
    assert hb.strand_capacity(hb.SID_MODE_HAMMING)==512
    assert list(book.values())[:512]==list(hb.codebook(512,mode).values())
    for sid in (512,1594,8191):
        assert rp.si_candidates(book[sid],mode,0)==((sid,8),)
        payload=('1011001'*15)[:103]
        encoded,*_=hb.select_encoding(payload,sid,book[sid])
        assert hb.decode_body(encoded.dna,sid,0,mode).payload_bits==payload


def test_new_rinf_is_reversible_and_preserves_old():
    data=b'\x00'*100+b'\xff'*100+bytes(range(256))
    old,_,_=lc.encode_bytes(data,'RINF_B16')
    dna,_,meta=lc.encode_bytes(data,'RINF_SHUFFLED_B16')
    assert dna!=old and len(dna)==len(old)
    assert lc.decode_dna(dna,'RINF_SHUFFLED_B16',meta)[0]==data
    assert lc.decode_dna(old,'RINF_B16',meta)[0]==data
    assert lc.encode_bytes(data,'RINF_SHUFFLED_B16')[0]==dna


def test_stage_progress_is_global_and_monotonic():
    strands,_=hb.encode_bytes(b'Two strands to check global stage ordering')
    updates=[]
    run=pm.simulate([s.row() for s in strands],pm.Settings(),stage_progress=lambda stage,p:updates.append((stage,p)))
    names=[stage for i,(stage,p) in enumerate(updates) if i==0 or stage!=updates[i-1][0]]
    assert names==['Synthesis','Storage','PCR','Sequencing']
    for name in names:
        values=[p for stage,p in updates if stage==name]
        assert values[0]==0 and values[-1]==1 and values==sorted(values)
    assert run['settings']['syn_sub']==.0001
