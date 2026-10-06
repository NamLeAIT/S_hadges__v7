import json
import subprocess
import sys
import random
from pathlib import Path
from dataclasses import asdict,replace
import zlib
import pytest
import hedges_bridge as hb
import process_model as pm
import read_pipeline as rp


def test_decode_fresh_process_without_original(tmp_path):
    data=b'DNA independent run'
    strands,_=hb.encode_bytes(data)
    sim=pm.simulate([s.row() for s in strands],pm.Settings(synthesis=False,storage=False,pcr=False,sequencing=False,coverage=3))
    (tmp_path/'reads.fastq').write_bytes(sim['r1'])
    config=rp.Settings(expected_strands=len(strands),stored_bytes=len(data))
    (tmp_path/'public.json').write_text(json.dumps(asdict(config)))
    script=Path(__file__).resolve().parents[1]/'decode_fastq.py'
    out=subprocess.run([sys.executable,str(script),'--r1',str(tmp_path/'reads.fastq'),'--settings',str(tmp_path/'public.json'),'--output',str(tmp_path/'decoded.bin')],cwd=tmp_path,capture_output=True,text=True,timeout=60)
    assert out.returncode==0,out.stderr
    assert (tmp_path/'decoded.bin').read_bytes()==data


def test_exported_record_order_and_names_irrelevant():
    data=b'order independent data'
    strands,_=hb.encode_bytes(data)
    sim=pm.simulate([s.row() for s in strands],pm.Settings(synthesis=False,storage=False,pcr=False,sequencing=False,coverage=4))
    lines=sim['r1'].decode().splitlines();units=[lines[i:i+4] for i in range(0,len(lines),4)]
    random.Random(77).shuffle(units)
    for i,u in enumerate(units):u[0]='@unrelated_'+str(i)
    changed=('\n'.join(x for u in units for x in u)+'\n').encode()
    cfg=rp.Settings(expected_strands=len(strands),stored_bytes=len(data))
    a=rp.run(sim['r1'],cfg);b=rp.run(changed,cfg)
    assert a['data']==b['data']==data
    assert a['consensus_rows']==b['consensus_rows']


def test_noisy_si_is_correct_or_rejected_not_header_rescued():
    data=b'DNA test';strands,_=hb.encode_bytes(data);s=strands[0].full_strand
    # SI deletion; public candidate search may accept or report ambiguity.
    noisy=s[:23]+s[24:]
    raw=('@strand_0\n'+noisy+'\n+\n'+'I'*len(noisy)+'\n').encode()
    out=rp.run(raw,rp.Settings(expected_strands=1,stored_bytes=8,joint_si_body=True))
    assert out['data'] in (None,data)


def test_quality_mode_controls_sequencing_not_inherited_flags():
    strands,_=hb.encode_bytes(b'DNA test')
    cfg=pm.Settings(synthesis=False,pcr=False,sequencing=True,quality_mode='Q-score',q_fixed=60,seq_sub=1,seq_ins=0,seq_del=0,coverage=10)
    sim=pm.simulate([s.row() for s in strands],cfg)
    # Error-rate slider is ignored in Q-score mode, as stated publicly.
    assert all(r['seq']==strands[0].full_strand for r in rp.parse(sim['r1']))


def test_lossless_text_wrapper_fix_retained():
    data=zlib.compress(b'compressed text example')
    strands,_=hb.encode_bytes(data)
    sim=pm.simulate([s.row() for s in strands],pm.Settings(synthesis=False,pcr=False,sequencing=False,coverage=3))
    out=rp.run(sim['r1'],rp.Settings(expected_strands=len(strands),stored_bytes=len(data)))
    assert out['data']==data
    assert out['restored_data']==b'compressed text example'
