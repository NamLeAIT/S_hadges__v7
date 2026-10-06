from pathlib import Path
import hashlib
from streamlit.testing.v1 import AppTest

APP=str(Path(__file__).resolve().parents[1]/'app.py')


def button(app,label):return next(x for x in app.button if x.label==label)


def test_single_page_fresh():
    app=AppTest.from_file(APP,default_timeout=120).run()
    assert not app.exception
    app.radio(key='consensus_source').set_value('Uploaded FASTQ').run()
    assert not app.exception
    assert not app.sidebar.radio
    text='\n'.join(x.value for x in app.markdown)
    assert all(label in text for label in ['1. Input','2. Compression','3. Encoding','4. NGS Strand Preparation','5. FASTQ Generation','6. Consensus','7. Decoding','8. Summarization'])
    assert all('dropout' not in str(x.label).lower() and 'reed' not in str(x.label).lower() for x in app.number_input)


def test_dt4dds_art_ui_through_process_and_decode():
    import pytest,art_sequencing,dt4dds_backend
    if not art_sequencing.detect_art_executable():pytest.skip('Real ART required')
    if not dt4dds_backend.available():pytest.skip('DT4DDS 1.1.0 required')
    app=AppTest.from_file(APP,default_timeout=120)
    # Inject the same source-file state the original upload panel would set.
    app.session_state['input_bytes']=b'DNA test'
    app.session_state['input_name']='test.txt'
    app.session_state['input_mime']='text/plain'
    app.session_state['input_sha']=hashlib.sha256(b'DNA test').hexdigest()
    app.run();assert not app.exception
    button(app,'Run Compression').click().run();assert not app.exception
    button(app,'Run DNA Encoding').click().run();assert not app.exception
    button(app,'Run Strand Preparation').click().run();assert not app.exception
    assert all(x.key!='channel_model' for x in app.radio)
    app.number_input(key='experimental_coverage').set_value(5).run()
    assert app.radio(key='experimental_sampling').value=='Pool distribution'
    app.radio(key='experimental_sampling').set_value('Equal coverage').run()
    button(app,'Run simulation').click().run();assert not app.exception
    run=app.session_state['process_run']
    assert run['engine']=='DT4DDS' and run['art']['engine']=='ART'
    assert run['settings']['sampling']=='Equal per strand'
    assert run['preset']['dt4dds_version']=='1.1.0'
    assert [x['Stage'] for x in run['stages']]==['Synthesis','Storage','PCR','Sampling','Sequencing']
    assert all(x['name'].startswith('read') and 'sid' not in x['name'].lower() for x in __import__('read_pipeline').parse(run['r1']))
    labels={x.label for x in app.metric}
    assert {'Sampled strands','Missing strands','Mean coverage','Error rate'}<=labels
    text='\n'.join(x.value for x in app.markdown)
    assert 'Simulation results' in text
    assert 'Scientific profile' not in text
    assert 'Accumulated error burden at each checkpoint' not in text
    assert 'Molecular abundance and physical sampling' not in text
    button(app,'Build consensus').click().run();assert not app.exception
    assert app.session_state['sim_decode']['result']['data'] is None
    button(app,'Decode consensus').click().run();assert not app.exception
    assert app.session_state['sim_decode']['result']['data']==b'DNA test'
    assert all('Completed · saved run' not in x.value for x in app.markdown)
    app.number_input(key='experimental_coverage').set_value(6).run();assert not app.exception
    assert 'process_run' not in app.session_state
    # Uploaded reads must not inherit the generated file's private layout.
    app.radio(key='consensus_source').set_value('Uploaded FASTQ').run()
    if 'upload_stored_bytes' in app.session_state:assert app.session_state['upload_stored_bytes']==0


def test_previous_step_consensus_decode_and_three_column_summary():
    import json
    from unittest.mock import patch
    import hedges_bridge as hb
    import process_model as pm
    raw=b'DNA test'
    strands,meta=hb.encode_bytes(raw,filename='test.txt',mime='text/plain')
    rows=[s.row() for s in strands]
    run=pm.simulate(rows,pm.Settings(synthesis=False,pcr=False,sequencing=False))
    run.update(design_signature=hashlib.sha256(json.dumps(rows,indent=2,ensure_ascii=False).encode()).hexdigest(),stored_reference=raw)
    app=AppTest.from_file(APP,default_timeout=120)
    for k,v in dict(input_bytes=raw,input_name='test.txt',input_mime='text/plain',stored_bytes=raw,
                    stored_name='test.txt',stored_mime='text/plain',dna_design='TVH-V5.2 HEDges',
                    encoded_ready=True,hedges_strands=strands,hedges_meta=meta,strand_rows=rows,
                    compression_method='No compression',process_run=run).items():app.session_state[k]=v
    # This UI test isolates source routing; real ART is tested separately.
    with patch('channel_ui.render',return_value=None):
        app.run();assert not app.exception
        button(app,'Build consensus').click().run();assert not app.exception
        button(app,'Decode consensus').click().run();assert not app.exception
        assert app.session_state['sim_decode']['result']['data']==raw
        text='\n'.join(x.value for x in app.markdown)
        assert all(x in text for x in ['Original','Compressed / encoded','Decoded','Encode–decode analysis'])
        assert any('Stored-byte integrity' in str(x.value) and 'PASS' in str(x.value) for x in app.dataframe)
        app.radio(key='consensus_source').set_value('Uploaded FASTQ').run();assert not app.exception
        assert not any('Stored-byte integrity' in str(x.value) for x in app.dataframe)


def test_import_package_then_consensus_then_decode():
    import io,json,zipfile
    from dataclasses import asdict
    from unittest.mock import patch
    import streamlit as st
    import hedges_bridge as hb
    import process_model as pm
    import read_pipeline as rp
    import pandas as pd
    strands,_=hb.encode_bytes(b'DNA test');rows=[s.row() for s in strands]
    run=pm.simulate(rows,pm.Settings(synthesis=False,pcr=False,sequencing=False))
    package=io.BytesIO()
    with zipfile.ZipFile(package,'w') as z:
        z.writestr('R1.fastq',run['r1'])
        z.writestr('decoder_settings.json',json.dumps(asdict(rp.Settings(expected_strands=1,stored_bytes=8))))
        z.writestr('evaluation/reference_strands.csv',pd.DataFrame(rows).to_csv(index=False))
    package.name='sequencing_output.zip'
    real=st.file_uploader
    def upload(label,*a,**k):
        if label=='FASTQ / R1 or sequencing ZIP':return package
        return real(label,*a,**k)
    with patch.object(st,'file_uploader',upload):
        app=AppTest.from_file(APP,default_timeout=120).run()
        app.radio(key='consensus_source').set_value('Uploaded FASTQ').run();assert not app.exception
        button(app,'Load public settings from ZIP').click().run();assert not app.exception
        button(app,'Build consensus').click().run();assert not app.exception
        assert app.session_state['upload_decode']['result']['data'] is None
        button(app,'Decode consensus').click().run();assert not app.exception
        assert app.session_state['upload_decode']['result']['data']==b'DNA test'


def test_single_workflow_and_missing_art_fails_closed():
    app=AppTest.from_file(APP,default_timeout=120)
    app.session_state['input_bytes']=b'DNA test'
    app.session_state['input_name']='test.txt';app.session_state['input_mime']='text/plain'
    app.session_state['input_sha']=hashlib.sha256(b'DNA test').hexdigest()
    app.run();button(app,'Run Compression').click().run()
    button(app,'Run DNA Encoding').click().run();button(app,'Run Strand Preparation').click().run()
    assert all(x.key!='channel_model' for x in app.radio)
    assert button(app,'Run simulation')
    app.text_input(key='experimental_art_executable').set_value('/missing/art').run()
    assert button(app,'Run simulation').disabled
    assert 'process_run' not in app.session_state
    assert any('not found' in x.value for x in app.info)
