from exports import collect
"""Single-page HEDGES workbench; FASTQ source and evaluation stay separate."""
import json
import hashlib
import io
import zipfile
import html
from pathlib import Path
from dataclasses import asdict, replace
import pandas as pd
import streamlit as st
import original_ui as ui
import hedges_bridge as hb
import process_model as pm
import read_pipeline as rp
import visual_ui as vu


def js(value):return json.dumps(value,indent=2,ensure_ascii=False).encode()
def digest(value):return hashlib.sha256(value).hexdigest()
def csv(rows):return pd.DataFrame(rows).to_csv(index=False).encode()
def zipped(files):
    b=io.BytesIO()
    with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as z:
        for name,value in files.items():z.writestr(name,value)
    return b.getvalue()


def current_public_settings(rows):
    codec=st.session_state['dna_design']
    return rp.Settings(codec=codec,sid_mode=(st.session_state.get('hedges_meta') or {}).get('sid_mode',hb.SID_MODE_HAMMING),
                       expected_strands=len(rows),stored_bytes=0 if codec=='TVH-V5.2 HEDges' else len(st.session_state.get('stored_bytes') or b''),dna_length=len(st.session_state.get('legacy_dna') or ''))


def settings_widgets(prefix,consensus_only=False):
    st.markdown('#### Consensus settings' if consensus_only else '#### Decoder settings')
    imported=st.file_uploader('Load public decoder settings (optional)',type='json',key=prefix+'_settings_file')
    field_keys=[k for k in rp.Settings.__dataclass_fields__ if k not in ('fbr','rbr')]
    if imported and st.button('Apply loaded settings',key=prefix+'_apply'):
        try:
            obj=asdict(rp.load_settings(imported.getvalue()))
            if obj['codec']!='TVH-V5.2 HEDges':raise ValueError('Only HEDGES files are supported.')
            for name in field_keys:st.session_state[prefix+'_'+name]=obj[name]
            st.success('Public settings loaded. Review them below.')
        except Exception as exc:st.error(str(exc))
    codec='TVH-V5.2 HEDges'
    st.write('**HEDGES · TVH-V5.2**')
    mode=st.radio('SI coding',[hb.SID_MODE_HAMMING,hb.SID_MODE_NO_HAMMING],horizontal=True,key=prefix+'_sid_mode')
    method=st.radio('Consensus method',['Quality-weighted voting','Majority voting'],horizontal=True,key=prefix+'_method')
    with st.expander('Alignment and quality settings'):
        minq=st.slider('Minimum mean Q-score',0,60,0,key=prefix+'_min_q')
        primer=st.slider('Allowed edits per primer',0,3,2,key=prefix+'_primer_edits')
        si=st.slider('Allowed SI edits (Hamming codebook)',0,2,1,key=prefix+'_si_edits')
        joint=st.checkbox('Check ambiguous SI boundaries with HEDGES',False,key=prefix+'_joint_si_body',help='At most 8 candidates. Accept only one verified boundary; otherwise report ambiguity. This is not an arbitrary-indel guarantee.')
        support=st.slider('Minimum consensus support',.51,1.,.6,.01,key=prefix+'_support')
        minimum=st.number_input('Minimum reads per strand',1,100,1,key=prefix+'_min_reads')
        maximum=st.number_input('Maximum read units used per strand',1,500,500,key=prefix+'_max_reads')
        outlier=st.slider('Maximum relative distance from cluster medoid',.05,1.,.25,.05,key=prefix+'_outlier_fraction')
        overlap=st.number_input('Minimum paired-end overlap (nt)',10,150,25,key=prefix+'_min_overlap')
        overlap_error=st.slider('Maximum overlap error fraction',0.,.5,.15,.01,key=prefix+'_overlap_error')
    if not consensus_only:
        with st.expander('HEDGES and file layout'):
            radius=st.slider('Maximum HEDGES edit radius (K)',0,4,1,key=prefix+'_max_edits')
            n=st.number_input('Expected strands (0 = unspecified)',0,65536,0,key=prefix+'_expected_strands')
            size=st.number_input('Stored file bytes (0 = infer framing for HEDGES)',0,100000000,0,key=prefix+'_stored_bytes')
            dna=0
        with st.expander('Candidate search · uncertain consensus'):
            rescue=st.checkbox('Try alternatives for uncertain positions',True,key=prefix+'_rescue_uncertain')
            uncertain=st.number_input('Maximum uncertain positions',0,12,6,key=prefix+'_max_uncertain')
            candidates=st.number_input('Maximum candidate combinations',1,4096,256,key=prefix+'_max_candidates')
            seconds=st.number_input('Time budget per strand (seconds)',.1,600.,30.,key=prefix+'_search_seconds')
            st.write('Incomplete search = SEARCH_LIMIT. No candidate is accepted until the declared search is complete.')
    else:
        defaults=rp.Settings()
        radius=st.session_state.get(prefix+'_max_edits',defaults.max_edits)
        n=st.session_state.get(prefix+'_expected_strands',defaults.expected_strands)
        size=st.session_state.get(prefix+'_stored_bytes',defaults.stored_bytes)
        dna=st.session_state.get(prefix+'_dna_length',defaults.dna_length)
        rescue=st.session_state.get(prefix+'_rescue_uncertain',defaults.rescue_uncertain)
        uncertain=st.session_state.get(prefix+'_max_uncertain',defaults.max_uncertain)
        candidates=st.session_state.get(prefix+'_max_candidates',defaults.max_candidates)
        seconds=st.session_state.get(prefix+'_search_seconds',defaults.search_seconds)
    cfg=rp.Settings(codec=codec,sid_mode=mode,method=method,min_q=minq,primer_edits=primer,si_edits=si,joint_si_body=joint,
                    support=support,min_reads=int(minimum),max_reads=int(maximum),outlier_fraction=outlier,min_overlap=int(overlap),
                    overlap_error=overlap_error,max_edits=radius,expected_strands=int(n),stored_bytes=int(size),dna_length=int(dna),
                    rescue_uncertain=rescue,max_uncertain=int(uncertain),max_candidates=int(candidates),search_seconds=float(seconds))
    collect('Download displayed decoder settings',js(asdict(cfg)),'decoder_settings.json','application/json',key=prefix+'_export')
    return cfg


def show_result(result,prefix,reference=None):
    if result.get('status')=='CONSENSUS_READY':
        st.info('Consensus is ready. Run the decoder below to recover the file.');return
    if result.get('data') is not None:
        st.success('File reconstructed. Check the hash below to establish exact equality.')
        st.write('Stored SHA-256:',result.get('sha256'))
        if reference is not None:
            if result['data']==reference:st.success('Stored-byte comparison: exact match (SHA-256 PASS).')
            else:st.error('Stored-byte comparison: mismatch. The reconstructed candidate is not exact.')
        data=result.get('restored_data',result['data'])
        name=result.get('restored_name') or result.get('filename') or 'recovered'+(result.get('extension') or '.bin')
        mime=result.get('restored_mime') or result.get('mime') or 'application/octet-stream'
        ui._preview_bytes(data,name,mime,'Reconstructed file',key=prefix+'_preview')
        collect('Download reconstructed file',data,name,mime,key=prefix+'_download')
    else:st.error(result.get('error') or f'Reconstruction: {result.get("status","FAIL")}')
    stats=result.get('read_stats',{})
    if stats:
        a,b,c=st.columns(3)
        a.metric('Read units',stats['input_units']);b.metric('Assigned units',stats['assigned_units']);c.metric('Consensus strands',len(result.get('consensus',{})))
    with st.expander('Consensus and decoding details',expanded=False):
        st.dataframe(pd.DataFrame(result.get('consensus_rows',[])),hide_index=True,width='stretch')
        st.dataframe(pd.DataFrame(result.get('strand_results',[])),hide_index=True,width='stretch')
        st.json(stats)
        collect('Download consensus positions',csv(result.get('positions',[])),'consensus_positions.csv',key=prefix+'_positions')
    report={k:v for k,v in result.items() if k not in ('data','restored_data','assignments','observed_bodies')}
    collect('Download result report',js(report),'decode_report.json',key=prefix+'_report')


def prep():
    with st.container(border=True):
        ui.step_header(4,'NGS Strand Preparation')
        if not st.session_state.get('encoded_ready'):
            st.info('Complete DNA Encoding first.');return []
        left,right=st.columns(2,gap='large')
        with left:
            st.markdown('#### Prepare clean strands')
            if st.button('Run Strand Preparation',type='primary',width='stretch'):
                rows=ui._ensure_strands();st.session_state['strand_rows']=rows
                for key in ('process_run','sim_decode'):st.session_state.pop(key,None)
            rows=st.session_state.get('strand_rows') or []
            if rows:
                a,b=st.columns(2);a.metric('Strands',len(rows));b.metric('Total nt',sum(len(r['Full strand']) for r in rows))
                collect('Download prepared strands CSV',ui._df_csv(rows),'prepared_strands.csv','text/csv')
                collect('Download full strands TXT',('\n'.join(r['Full strand'] for r in rows)+'\n').encode(),'full_strands.txt')
                settings=current_public_settings(rows)
                collect('Download public encoding settings',js(asdict(settings)),'encoding_settings.json','application/json')
        with right:
            st.markdown('#### Strand structure')
            if rows:
                sid=st.number_input('Inspect clean strand',0,len(rows)-1,0)
                ui.show_segmented(rows[sid])
                st.dataframe(pd.DataFrame([{k:r.get(k) for k in ('Strand ID','Payload length','Full strand length','GC content','Longest homopolymer')} for r in rows]),hide_index=True,width='stretch')
            else:st.info('Prepare strands to see their structure.')
        return rows


def processes(rows):
    import channel_ui
    channel_ui.render(rows,ui,digest,js,csv,zipped,current_public_settings)


def decode_panel(prefix='upload'):
    import imports as inputs
    st.markdown('### 1 · Import sequencing reads')
    f1=st.file_uploader('Sequencing package or FASTQ / R1',type=['zip','fastq','fq','gz'],key='uploaded_r1')
    f2=st.file_uploader('R2 — paired-end only',type=['fastq','fq','gz'],key='uploaded_r2')
    files={};r1=None;r2=b'';rows=None;expected=''
    if f1:
        try:
            if f1.name.lower().endswith('.zip'):
                files=inputs.package(f1.getvalue());r1=files['R1.fastq'];r2=files.get('R2.fastq',b'')
                if f2:st.warning('R2 comes from the ZIP. The separate R2 upload is ignored.')
            else:r1=f1.getvalue();r2=f2.getvalue() if f2 else b''
        except Exception as exc:st.error(str(exc));return
    if 'decoder_settings.json' in files:
        if st.button('Apply public settings from package',type='primary'):
            try:
                public=rp.load_settings(files['decoder_settings.json'])
                for k,v in asdict(public).items():st.session_state[prefix+'_'+k]=v
                st.success('Public settings applied. Review in Advanced settings.')
            except Exception as exc:st.error(str(exc))
    with st.expander('Advanced settings · consensus and HEDGES'):
        cfg=settings_widgets(prefix)
    with st.expander('Reference comparison · evaluation only',expanded='evaluation/reference_strands.csv' in files):
        rf=st.file_uploader('Reference strands CSV (optional)',type='csv',key='evaluation_reference')
        refraw=rf.getvalue() if rf else files.get('evaluation/reference_strands.csv')
        compare=st.checkbox('Compare results with reference after processing',value=True,key='evaluate_reference') if refraw else False
        if compare and refraw:
            try:rows=inputs.references(refraw)
            except Exception as exc:st.error(str(exc))
        expected=st.text_input('Expected stored SHA-256 (optional)',value=files.get('evaluation/expected_hash.txt',b'').decode().strip(),key='upload_hash_check')
    if r1:
        try:
            raw=rp.parse(r1)+(rp.parse(r2) if r2 else [])
            vu.length_summary([len(r['seq']) for r in raw],'Imported FASTQ lengths (nt)')
            vu.fastq_preview(r1,r2)
        except Exception as exc:st.error(str(exc));return
    process_reads(r1,r2,cfg,prefix,rows,expected)


def process_reads(r1,r2,cfg,prefix,rows=None,expected=''):
    import candidate_search
    left,right=st.columns([1,1.35],gap='large')
    with left:
        if cfg is None:
            with st.expander('Public design and alignment',expanded=True):
                cfg=settings_widgets(prefix,consensus_only=True)
        consensus_cfg=asdict(cfg)
        for key in ('rescue_uncertain','max_uncertain','max_candidates','search_seconds'):consensus_cfg.pop(key)
        if not cfg.joint_si_body:consensus_cfg.pop('max_edits')
        fingerprint=digest((r1 or b'')+r2+js(consensus_cfg))
        if st.button('Build consensus',type='primary',disabled=not r1,key=prefix+'_build'):
            try:
                bar=st.progress(0.)
                result=rp.run(r1,cfg,r2,bar.progress,decode_file=False)
                st.session_state[prefix+'_decode']={'fingerprint':fingerprint,'result':result}
            except Exception as exc:st.session_state.pop(prefix+'_decode',None);st.error(str(exc))
    saved=st.session_state.get(prefix+'_decode')
    if not saved:
        with right:st.info('Build consensus to compare similarity before and after.')
        ui.step_header(7,'Decoding');st.info('Build consensus first.')
        ui.step_header(8,'Summarization');st.info('Results appear after decoding.')
        return
    if saved['fingerprint']!=fingerprint:
        st.warning('Input or settings changed. Build consensus again.')
        ui.step_header(7,'Decoding');st.info('Rebuild consensus first.')
        ui.step_header(8,'Summarization');st.info('Results appear after decoding.')
        return
    result=saved['result']
    with right:
        vu.consensus_summary(result,rows)
        capped=sum(r.get('capped',0) for r in result.get('consensus_rows',[]))
        if capped:st.warning(f'{capped:,} assigned read units exceeded the consensus cap and were not used in voting. Increase the cap in Alignment and quality settings (maximum 500).')
    with st.expander('Consensus strands, errors and lengths',expanded=False):
        vu.consensus_results(result,prefix,rows)
        vu.length_comparison(result)
    collect('Consensus positions',csv(result.get('positions',[])),'consensus_positions.csv')
    collect('Consensus sequences',''.join(f'>assigned_strand_{sid}\n{seq}\n' for sid,seq in result['all_consensus'].items()).encode(),'consensus.fasta')
    ui.step_header(7,'Decoding')
    with st.expander('Advanced · HEDGES search and file layout'):
        radius=st.slider('Maximum HEDGES edit radius (K)',0,4,cfg.max_edits,key=prefix+'_max_edits')
        rescue=st.checkbox('Try alternatives for uncertain positions',cfg.rescue_uncertain,key=prefix+'_rescue_uncertain')
        uncertain=st.number_input('Maximum uncertain positions',0,12,cfg.max_uncertain,key=prefix+'_max_uncertain')
        candidates=st.number_input('Maximum candidate combinations',1,4096,cfg.max_candidates,key=prefix+'_max_candidates')
        seconds=st.number_input('Time budget per strand (seconds)',.1,600.,cfg.search_seconds,key=prefix+'_search_seconds')
        n=st.number_input('Expected strands',0,65536,cfg.expected_strands,key=prefix+'_expected_strands')
        size=st.number_input('Stored file bytes',0,100000000,cfg.stored_bytes,key=prefix+'_stored_bytes')
    cfg=replace(cfg,max_edits=radius,rescue_uncertain=rescue,max_uncertain=int(uncertain),max_candidates=int(candidates),search_seconds=float(seconds),expected_strands=int(n),stored_bytes=int(size),dna_length=0)
    decode=st.button('Decode consensus',type='primary',key=prefix+'_decode_button')
    retry=st.button('Retry limited strands',key=prefix+'_retry',disabled=not any(r.get('status')=='SEARCH_LIMIT' for r in result.get('search_results',[])))
    if decode or retry:
        with st.spinner('Searching candidates and decoding…'):
            bar=st.progress(0.)
            output=candidate_search.decode(result,cfg,bar.progress,retry_only=retry)
            for k in ('data','restored_data','restored_name','restored_mime','sha256','error'):result.pop(k,None)
            result.update(output)
            result['settings']=asdict(cfg)
    if result.get('search_results'):
        st.dataframe(pd.DataFrame(result['search_results']).drop(columns=['payload_bits'],errors='ignore'),hide_index=True,width='stretch')
    if result['status']!='CONSENSUS_READY':
        show_result(result,prefix)
        if expected:
            if result.get('sha256')==expected.strip().lower():st.success('Exact stored file · SHA-256 PASS')
            else:st.error('Stored SHA-256 mismatch or file not reconstructed')
        if rows and result.get('strand_results'):
            measured=vu.ev.hedges_report(rows,result)
            st.dataframe(pd.DataFrame(measured),hide_index=True,width='stretch')
            collect('HEDGES evaluation',csv(measured),'hedges_evaluation.csv')
    summarization(result,prefix,rows)
    collect('Decoder report',js({k:v for k,v in result.items() if k not in ('data','restored_data')}),'decode_report.json')


def download_center(page):
    st.markdown('### 📁 All results')
    sections=st.session_state.setdefault('export_sections',{})
    sections['simulation' if page.startswith('1') else 'consensus']=dict(st.session_state.get('all_exports',{}))
    files={section+'/'+name:value for section,items in sections.items() for name,value in items.items()}
    if files:
        with st.expander('Reports, DNA sequences and recovered files'):
            st.write('Reports and recovered files')
            st.write(sorted(files))
            st.download_button('Download All results',zipped({'all/'+name:value for name,value in files.items()}),'All_results.zip')
    else:st.write('Results will appear here after processing.')


def summarization(result,prefix,rows):
    """Reference comparisons are post-processing only, never decoder inputs."""
    with st.container(border=True):
        ui.step_header(8,'Summarization')
        st.markdown('#### File comparison')
        local=prefix=='sim'
        original=st.session_state.get('input_bytes') if local else None
        stored=st.session_state.get('stored_bytes') if local else None
        candidates=st.session_state.get('compression_rows',[]) if local else []
        selected=next((c for c in candidates if c['Method']==st.session_state.get('compression_method')),None)
        restored=result.get('restored_data')
        if restored is None:restored=result.get('data')
        previews=[('Original',original,st.session_state.get('input_name','input.bin'),st.session_state.get('input_mime')),
                  ('Compressed / encoded',stored,st.session_state.get('stored_name','stored.bin'),st.session_state.get('stored_mime')),
                  ('Decoded',restored,result.get('restored_name') or result.get('filename') or 'recovered.bin',result.get('restored_mime') or result.get('mime'))]
        for i,(col,item) in enumerate(zip(st.columns(3,gap='large'),previews)):
            label,data,name,mime=item
            with col:
                st.markdown('##### '+label)
                if data is None:
                    st.info('Not available' if i<2 else 'No reconstructed file')
                    continue
                preview_data,preview_name,preview_mime=data,name,mime
                if i==1 and selected:preview_data,preview_name,preview_mime=ui._candidate_preview(selected)
                ui._preview_bytes(preview_data,preview_name,preview_mime,'Preview',key=prefix+'_summary_'+str(i))
                props=[{'Property':'Extension','Value':Path(preview_name).suffix or '—'},
                       {'Property':'Type','Value':preview_mime or 'application/octet-stream'},
                       {'Property':'Stored size' if i==1 else 'Size','Value':f'{len(data)/1024:.2f} KB'}]
                if preview_mime and preview_mime.startswith('image/') and ui.Image:
                    try:
                        with ui.Image.open(io.BytesIO(preview_data)) as im:props.append({'Property':'Image size','Value':f'{im.width} × {im.height} px'})
                    except Exception:pass
                ui.property_table(props)
        st.markdown('#### Compression analysis')
        analysis=st.session_state.get('compression_analysis_rows',[]) if local else []
        chosen=next((r for r in analysis if r['Method']==st.session_state.get('compression_method')),None)
        if chosen:
            row={'Method':ui._format_method_name(chosen['Method']),'Original data':f'{len(original):,} B',
                 'Original extension':Path(st.session_state.get('input_name','')).suffix,
                 'Compressed data':f'{len(stored):,} B','Compressed extension':chosen.get('Format'),
                 'Compression ratio':f"{chosen['Compression ratio']:.2f}×",'Size reduction':f"{chosen['Size reduction (%)']:.2f}%"}
            for key in ('PSNR','SSIM','MAE','Audio SNR','Audio correlation','Audio NRMSE','Video SSIM','Video PSNR'):
                if key in chosen:row[key]=chosen[key]
            st.dataframe(pd.DataFrame([row]),hide_index=True,width='stretch')
        else:st.info('Compression details are available for files encoded in this session.')
        st.markdown('#### Encode–decode analysis')
        strands=rows or []
        bodies=[str(r.get('Payload','')) for r in strands]
        full=[str(r.get('Full strand','')) for r in strands]
        nt=sum(map(len,bodies));fullnt=sum(map(len,full))
        properties=[{'Property':'DNA design','Value':'HEDGES · TVH-V5.2'},
                    {'Property':'Strand format','Value':'FBR + SI + HEDGES BODY + RBR'},
                    {'Property':'BODY DNA length','Value':f'{nt:,} nt' if strands else '—'},
                    {'Property':'BODY GC content','Value':f'{100*sum(s.count("G")+s.count("C") for s in bodies)/max(1,nt):.2f}%' if strands else '—'},
                    {'Property':'BODY longest homopolymer','Value':max([ui.legacy_codec.longest_homopolymer(s) for s in bodies] or [0]) if strands else '—'},
                    {'Property':'Designed strands','Value':len(strands) if strands else result['settings'].get('expected_strands') or 'Unspecified'},
                    {'Property':'Full strand DNA length','Value':f'{fullnt:,} nt' if strands else '—'},
                    {'Property':'Mean strand GC content','Value':f'{100*sum(s.count("G")+s.count("C") for s in full)/max(1,fullnt):.2f}%' if strands else '—'},
                    {'Property':'Full strand longest homopolymer','Value':max([ui.legacy_codec.longest_homopolymer(s) for s in full] or [0]) if strands else '—'}]
        checked=vu.ev.hedges_report(strands,result) if strands and result.get('strand_results') else []
        if checked:
            exact=sum(r['Payload comparison']=='EXACT' for r in checked)
            properties.append({'Property':'Exact payloads / designed strands','Value':f'{exact} / {len(strands)}'})
        recovered_stored=result.get('data')
        match=(recovered_stored==stored) if local and recovered_stored is not None else None
        properties.extend([{'Property':'Stored-byte integrity','Value':('PASS' if match else 'FAIL') if match is not None else ('Framing validated' if recovered_stored is not None else 'Not reconstructed')},
                           {'Property':'Recovery status','Value':result.get('status','—')}])
        ui.property_table(properties)
        a,b,c=st.columns(3)
        a.metric('Read units',result['read_stats']['input_units'])
        b.metric('Assigned units',result['read_stats']['assigned_units'])
        c.metric('Consensus strands',sum(bool(s) for s in result['all_consensus'].values()))
        if local and st.session_state.get('compression_lossy'):
            st.info('Recovery checks the selected compressed file, not equality with the original file before lossy compression.')
        collect('Summary',csv(properties),'summary.csv')


def consensus_sources(rows,run):
    import imports as inputs
    with st.container(border=True):
        ui.step_header(6,'Consensus')
        left,right=st.columns(2,gap='large')
        with left:
            st.markdown('#### Import from previous step')
            available=bool(run and rows and run['design_signature']==digest(js(rows)))
            if available:
                st.success('Generated FASTQ is ready')
                st.write(f"{len(rp.parse(run['r1'])):,} read units · {'Paired-end' if run['r2'] else 'Single-end'}")
            else:st.info('Generate FASTQ in Step 5, or upload on the right.')
        with right:
            st.markdown('#### Upload FASTQ')
            f1=st.file_uploader('FASTQ / R1 or sequencing ZIP',type=['zip','fastq','fq','gz'],key='uploaded_r1')
            f2=st.file_uploader('R2 (optional)',type=['fastq','fq','gz'],key='uploaded_r2')
        source=st.radio('Process reads from',['Previous step','Uploaded FASTQ'],horizontal=True,key='consensus_source')
        if source=='Previous step':
            prefix='sim';r1=run['r1'] if available else None;r2=run['r2'] if available else b''
            reference=rows if available else None;expected=digest(run['stored_reference']) if available else ''
            if available:
                public=current_public_settings(rows);signature=digest(js(asdict(public)))
                if st.session_state.get('sim_public_signature')!=signature:
                    for k,v in asdict(public).items():st.session_state['sim_'+k]=v
                    st.session_state['sim_public_signature']=signature
        else:
            prefix='upload';r1=None;r2=b'';reference=None;expected='';files={}
            if f1:
                try:
                    if f1.name.lower().endswith('.zip'):
                        files=inputs.package(f1.getvalue());r1=files['R1.fastq'];r2=files.get('R2.fastq',b'')
                    else:r1=f1.getvalue();r2=f2.getvalue() if f2 else b''
                except Exception as exc:st.error(str(exc))
            if 'decoder_settings.json' in files and st.button('Load public settings from ZIP'):
                try:
                    public=rp.load_settings(files['decoder_settings.json'])
                    if public.codec!='TVH-V5.2 HEDges':raise ValueError('Only HEDGES files are supported in this app.')
                    for k,v in asdict(public).items():st.session_state[prefix+'_'+k]=v
                except Exception as exc:st.error(str(exc))
            with st.expander('Reference comparison (optional)'):
                rf=st.file_uploader('Designed strands CSV · evaluation only',type='csv',key='evaluation_reference')
                raw=rf.getvalue() if rf else files.get('evaluation/reference_strands.csv')
                if raw:
                    try:reference=inputs.references(raw)
                    except Exception as exc:st.error(str(exc))
                expected=st.text_input('Expected stored SHA-256',value=files.get('evaluation/expected_hash.txt',b'').decode().strip(),key='upload_hash_check')
    process_reads(r1,r2,None,prefix,reference,expected)


def workflow_header(box):
    result=(st.session_state.get(('sim' if st.session_state.get('consensus_source','Previous step')=='Previous step' else 'upload')+'_decode') or {}).get('result',{})
    flags=[bool(st.session_state.get('input_bytes')),bool(st.session_state.get('stored_bytes')),
           bool(st.session_state.get('encoded_ready')),bool(st.session_state.get('strand_rows')),
           bool(st.session_state.get('process_run')),bool(result),result.get('data') is not None,result.get('data') is not None]
    labels=['Input','Compression','Encoding','NGS Strand Preparation','FASTQ Generation','Consensus','Decoding','Summarization']
    items=[]
    for n,(label,done) in enumerate(zip(labels,flags),1):
        state='Done' if done else 'Ready' if n==1 or flags[n-2] else 'Waiting'
        items.append(f'<div class="pipeline-step {"done" if done else ""}"><span class="step-num">{n}</span><span class="step-name">{html.escape(label)}</span><div class="step-state">{state}</div></div>')
    box.markdown('<style>.pipeline-steps{grid-template-columns:repeat(4,1fr)}@media(max-width:700px){.pipeline-steps{grid-template-columns:repeat(2,1fr)}}</style><div class="pipeline-steps">'+''.join(items)+'</div>',unsafe_allow_html=True)


def main():
    st.set_page_config(page_title=ui.APP_TITLE,page_icon='🧬',layout='wide')
    ui.apply_style()
    st.session_state['all_exports']={}
    with st.sidebar.expander('Quick guide'):
        st.write('Molecule: a physical copy of a strand. PCR cycle: one amplification round. R1/R2: reads from opposite ends of the same molecule. Quality: FASTQ line 4, Phred+33.')
        st.write('Generate FASTQ or upload your sequencing reads in Step 6. Only HEDGES encoding is supported.')
    st.markdown('<div class="hero-card"><div class="hero-title">🧬 DNA Data Storage System</div></div>',unsafe_allow_html=True)
    header=st.empty()
    ui.render_panel_1_upload();ui.render_panel_2_compression();ui.render_panel_3_encoding()
    rows=prep();processes(rows)
    consensus_sources(rows,st.session_state.get('process_run'))
    workflow_header(header)
    download_center('1. Workflow')


if __name__=='__main__':main()
