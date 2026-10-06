from dataclasses import replace
import evaluation as ev
import hedges_bridge as hb
import process_model as pm
import read_pipeline as rp


def test_edit_types_and_unread_ends():
    assert ev.edits('ACGT','ATGT')['S']==1
    assert ev.edits('ACGT','ACAGT')['I']==1
    assert ev.edits('ACGT','ACT')['D']==1
    assert ev.edits('ACGT','ACNT')['N']==1
    assert ev.edits('TTACGTAACC','ACGT',partial=True)['errors']==0
    assert ev.edits('TTACGTAACC','ACGT',partial=True)['reference_nt']==4


def test_no_false_reduction_from_fewer_sequences():
    rows=[{'Strand ID':0,'Payload':'ACGT'}]
    result={'consensus':{0:'ATGT'},'all_consensus':{0:'ATGT'},
            'observed_bodies':[{'sid':0,'sequence':'ATGT'} for _ in range(10)]}
    report=ev.consensus_report(rows,result)
    assert report['before_rate']==report['after_rate']==.25
    assert report['reduction_percent']==0 # 10 raw errors -> 1 consensus error is NOT 90% improvement
    result['all_consensus'][0]='ACGT'
    assert ev.consensus_report(rows,result)['reduction_percent']==100
    result['all_consensus'][0]=''
    missing=ev.consensus_report(rows,result)
    assert missing['matched_strands']==0 and missing['after_rate'] is None


def test_two_step_decode_and_evaluator_cannot_change_output():
    data=b'DNA test';strands,_=hb.encode_bytes(data);rows=[s.row() for s in strands]
    cfg=replace(pm.Settings(),synthesis=False,pcr=False,sequencing=False,coverage=3,paired=True)
    sim=pm.simulate(rows,cfg)
    counts=ev.sequencing(rows,sim)
    assert counts['totals']['errors']==0
    assert counts['coverage'][0]['Read units']==3
    assert counts['coverage'][0]['FASTQ records']==6
    c=rp.Settings(stored_bytes=8,expected_strands=1)
    result=rp.run(sim['r1'],c,sim['r2'],decode_file=False)
    assert result['status']=='CONSENSUS_READY' and result['data'] is None
    ev.consensus_report(rows,result)
    result.update(rp.reconstruct(result['consensus'],c))
    assert result['data']==data
    assert ev.hedges_report(rows,result)[0]['Payload comparison']=='EXACT'
    wrongrows=[{**rows[0],'Payload':'A'*76}]
    assert ev.consensus_report(wrongrows,result)['after_rate']>0
    assert result['data']==data
