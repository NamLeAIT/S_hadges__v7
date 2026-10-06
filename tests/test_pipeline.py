from dataclasses import replace,asdict
import json
import random
import hashlib
import pytest
import hedges_bridge as hb
import legacy_codec as lc
import process_model as pm
import read_pipeline as rp


def fq(rows):
    return ''.join(f'@arbitrary_{i}\n{r["Full strand"]}\n+\n'+('I'*len(r['Full strand']))+'\n' for i,r in enumerate(rows)).encode()


@pytest.fixture(scope='module')
def hedges():
    rows,meta=hb.encode_bytes(b'DNA test')
    return [r.row() for r in rows]


def clean(**kwargs):
    return replace(pm.Settings(synthesis=False,storage=False,pcr=False,sequencing=False,coverage=3),**kwargs)


def public(**kwargs):
    return replace(rp.Settings(expected_strands=1,stored_bytes=8),**kwargs)


def test_unchanged_hedges_clean(hedges):
    result=rp.run(fq(hedges),public())
    assert result['data']==b'DNA test'
    assert len(hedges[0]['Full strand'])==125


def test_synthesis_inheritance_fixed_coverage(hedges):
    cfg=clean(synthesis=True,syn_sub=.02,syn_ins=.005,syn_del=.005,pcr=True,pcr_sub=0,pcr_ins=0,pcr_del=0,efficiency=1,cycles=3,coverage=30)
    a=pm.simulate(hedges,cfg);b=pm.simulate(hedges,cfg)
    assert a==b
    assert len(rp.parse(a['r1']))==30
    assert all(t['sid']==0 for t in a['truth'])
    by_founder={}
    for t in a['truth']:by_founder.setdefault(t['founder'],set()).add(t['template'])
    assert all(len(v)==1 for v in by_founder.values())
    counts={s['Stage']:s for s in a['stages']}
    assert counts['PCR cycle 3']['molecules']==counts['Synthesis']['molecules']*8


def test_pcr_new_variants_without_dropout(hedges):
    sim=pm.simulate(hedges,clean(pcr=True,pcr_sub=.005,pcr_ins=.001,pcr_del=.001,cycles=4,efficiency=1))
    counts={s['Stage']:s for s in sim['stages']}
    assert counts['PCR cycle 4']['molecules']==5*16
    assert counts['PCR cycle 4']['variants']>5
    assert len(rp.parse(sim['r1']))==3


def test_storage_on_synthesis_output(hedges):
    sim=pm.simulate(hedges,clean(synthesis=True,syn_sub=0,syn_ins=0,syn_del=0,storage=True,storage_sub=0,storage_ct=1))
    assert all('C' not in t['sequence'] for t in sim['trace'] if t['stage']=='Storage')
    assert len(rp.parse(sim['r1']))==3


def test_fastq_no_truth_and_q_not_oracle(hedges):
    sim=pm.simulate(hedges,clean(synthesis=True,syn_sub=1,syn_ins=0,syn_del=0))
    reads=rp.parse(sim['r1'])
    assert all('strand' not in r['name'] and 'copy' not in r['name'] for r in reads)
    assert reads[0]['seq']!=hedges[0]['Full strand']
    assert set(reads[0]['q'])=={60} # inherited wrong bases can be Q60


def test_primer_substitution_and_indel(hedges):
    row=dict(hedges[0]);s=row['Full strand']
    for changed in (s[:6]+s[7:], s[:6]+'A'+s[6:], s[:6]+('C' if s[6]!='C' else 'G')+s[7:]):
        row['Full strand']=changed
        out=rp.run(fq([row]),public(primer_edits=2))
        assert out['data']==b'DNA test'


def test_reverse_orientation_and_adapter_flanks(hedges):
    row=dict(hedges[0]);row['Full strand']='AGGA'+pm.rc(row['Full strand'])+'TTGGAC'
    assert rp.run(fq([row]),public())['data']==b'DNA test'


def test_msa_corrects_indels_and_substitutions():
    s='ACGTTAGCACCTGACTAGCT';reads=[]
    for i in range(7):
        changed=s
        if i<3:changed=s[:4+i]+s[5+i:]
        elif i<6:changed=s[:10+i]+'A'+s[10+i:]
        reads.append({'seq':changed,'q':[35]*len(changed)})
    seq,positions,info=rp.consensus(reads,rp.Settings())
    assert seq==s
    assert info['uncertain']==0


def test_ties_are_unknown_not_a_guess():
    reads=[{'seq':'ACGT','q':[30]*4},{'seq':'AGGT','q':[30]*4}]
    s,p,info=rp.consensus(reads,rp.Settings(outlier_fraction=1))
    assert s=='ANGT' and info['uncertain']==1


def test_consensus_order_independent():
    reads=[{'seq':s,'q':[30]*len(s)} for s in ['ACGTTAC','ACGTAC','ACGTTAC','ACCTTAC']]
    a=rp.consensus(reads,rp.Settings())
    random.Random(32).shuffle(reads)
    assert rp.consensus(reads,rp.Settings())==a


def test_headers_not_used(hedges):
    raw=fq(hedges)
    assert rp.run(raw.replace(b'@arbitrary_0',b'@strand_400'),public())['data']==b'DNA test'
    fake=b'@strand_0\n'+b'A'*76+b'\n+\n'+b'I'*76+b'\n'
    out=rp.run(fake,public())
    assert out['data'] is None and out['read_stats']['assigned_units']==0


def test_pair_handling(hedges):
    sim=pm.simulate(hedges,clean(paired=True))
    out=rp.run(sim['r1'],public(),sim['r2'])
    assert out['data']==b'DNA test'
    assert out['read_stats']['assigned_units']==3 and out['read_stats']['raw_reads']==6


def test_short_overlap_pair(hedges):
    sim=pm.simulate(hedges,clean(paired=True,read_length=90))
    out=rp.run(sim['r1'],public(),sim['r2'])
    assert out['data']==b'DNA test'
    assert out['read_stats']['PAIRS_MERGED']==3


def test_bad_fastq_and_settings_rejected():
    with pytest.raises(ValueError):rp.parse(b'@r\nAC\n+\nI\n')
    with pytest.raises(ValueError):rp.load_settings(json.dumps({'clean_payload':'ACGT'}))
    with pytest.raises(ValueError):rp.validate(replace(rp.Settings(),support=.5))


def test_missing_strand_rejected(hedges):
    assert rp.run(fq(hedges),public(expected_strands=2))['data'] is None


def test_legacy_public_metadata_only():
    for codec in ['Simple Mapping','R0_B9','R1_B12','R2_B15','RINF_B16','RINF_SHUFFLED_B16']:
        data=b'Hello DNA!';dna,bits,meta=lc.encode_bytes(data,codec)
        rows=[]
        for sid,start in enumerate(range(0,len(dna),76)):
            si=''.join('ACGT'[(sid >> (2*(7-i)))&3] for i in range(8))
            body=dna[start:start+76];fill=('ACGT'*19)[:76-len(body)]
            rows.append({'Full strand':hb.FBR+si+body+fill+hb.RBR})
        c=rp.Settings(codec=codec,expected_strands=len(rows),stored_bytes=len(data),dna_length=len(dna))
        assert rp.run(fq(rows),c)['data']==data
