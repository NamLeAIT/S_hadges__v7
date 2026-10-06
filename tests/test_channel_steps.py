from dataclasses import replace
import process_model as pm
import evaluation as ev


def rows():
    return [{'Strand ID':0,'Payload':'ACGT'*10,'Full strand':'AAAA'+'ACGT'*10+'TTTT',
             'FBR':'AA','SI':'AA','RBR':'TTTT','Filler':''}]


def test_sequential_stages_equal_one_shot_and_do_not_mutate_checkpoint():
    import copy
    cfg=pm.Settings(founders=4,cycles=3,coverage=10,syn_sub=.05,storage=True,storage_sub=.03,
                    pcr_sub=.04,seq_sub=.02,paired=True)
    all_at_once=pm.simulate(rows(),cfg)
    checkpoint=None
    for name in ['Synthesis','Storage','PCR','Sequencing']:
        before=copy.deepcopy(checkpoint)
        new=pm.simulate(rows(),cfg,stop_after=name,checkpoint=checkpoint)
        assert checkpoint==before
        if name!='Sequencing':assert 'r1' not in new
        checkpoint=new
    for key in ['r1','r2','truth','stages','abundance','stage_quality','audit']:
        assert all_at_once[key]==checkpoint[key],key


def test_downstream_setting_does_not_change_saved_upstream_pool():
    cfg=pm.Settings(founders=2,cycles=1)
    syn=pm.simulate(rows(),cfg,stop_after='Synthesis')
    store=pm.simulate(rows(),replace(cfg,storage=False),stop_after='Storage',checkpoint=syn)
    assert syn['_populations']==store['_populations']
    pcr=pm.simulate(rows(),replace(cfg,efficiency=1.,pcr_sub=0.,pcr_ins=0.,pcr_del=0.),stop_after='PCR',checkpoint=store)
    assert sum(pcr['_populations'][0].values())==2*sum(store['_populations'][0].values())


def test_strict_similarity_and_read_span():
    assert ev.strict_similarity('ACGT','ACGA')==75.
    assert ev.strict_similarity('ACGT','ACG')==0.
    cfg=pm.Settings(synthesis=False,storage=False,pcr=False,sequencing=False,paired=True,read_length=20)
    run=pm.simulate(rows(),cfg)
    report=ev.sequencing(rows(),run)
    assert all(r['Strict similarity (%)']==100. for r in report['reads'])
