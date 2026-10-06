"""Stage-by-stage channel controls and evaluation; no truth enters decoding."""
from dataclasses import asdict, replace
import hashlib
import json
import html
import pandas as pd
import streamlit as st
import altair as alt
import process_model as pm
import hedges_bridge as hb
import evaluation as ev
import visual_ui as vu
from exports import collect

STAGES=['Synthesis','Storage','PCR','Sampling','Sequencing']
FIELDS=[['seed','scope','founders','synthesis','syn_sub','syn_ins','syn_del'],
        ['storage','storage_sub','storage_ct'],['pcr','cycles','efficiency','bias_sd','pcr_sub','pcr_ins','pcr_del','max_variants'],
        ['sequencing','coverage','sampling','paired','read_length','quality_mode','seq_sub','seq_ins','seq_del','q_profile','q_fixed','q_start','q_end','q_std']]

def event_table(run):
    return [{'Process':name,**{k:sum(r.get(k,0) for r in run['stages'] if r['Stage']==name or name=='PCR' and r['Stage'].startswith('PCR cycle')) for k in ['sub','ins','del','input_nt']}} for name in STAGES]

def chart_events(table,x='Process',order=None,title='Events / observed edits'):
    frame=pd.DataFrame(table).melt(id_vars=[x],value_vars=['sub','ins','del'],var_name='Error',value_name='Events')
    frame['Error']=frame['Error'].map({'sub':'Substitution','ins':'Insertion','del':'Deletion'})
    st.altair_chart(alt.Chart(frame).mark_bar().encode(x=alt.X(x+':N',sort=order or STAGES),y=alt.Y('Events:Q',title=title),color=alt.Color('Error:N',sort=['Substitution','Insertion','Deletion']),tooltip=[x,'Error','Events']),width="stretch")

def pipeline(box,states):
    icons=['🧬','🧊','🧪','💧','🔬'];verbs=['Twist synthesis','Storage exposure','Q5 amplification','Physical aliquot','ART HS25 reads']
    with box.container():
        cols=st.columns(4)
        for i,c in enumerate(cols):
            with c:
                with st.container(border=True):
                    st.markdown(f'**{icons[i]} 5.{i+1} {STAGES[i]}**')
                    st.write(verbs[i]+(' →' if i<3 else ''));st.write(states[i])

def before_after(run,name):
    samples=[r for r in run.get('audit',[]) if r['Stage']==name]
    st.markdown('**Before / after**')
    if not samples:
        st.info('No new mutation examples in this stage. Parents and unchanged copies remain in the pool.');return
    ids=sorted({r['Strand ID'] for r in samples})
    sid=st.selectbox('Strand ID',ids,key='audit_sid_'+name)
    group=[r for r in samples if r['Strand ID']==sid]
    idx=st.selectbox('Molecule example',range(len(group)),format_func=lambda i:f"{i+1} · founder {group[i]['Founder']}"+(f" · cycle {group[i]['Cycle']}" if group[i]['Cycle'] else ''),key='audit_example_'+name)
    r=group[idx];a,b,c=st.columns(3)
    a.metric('Before length (nt)',r['Before nt']);b.metric('After length (nt)',r['After nt']);c.metric('New events in this example',sum(r[k] for k in ['New S','New I','New D']))
    aa,bb=ev.alignment(r['Before'],r['After']);pos=0;out=[]
    for x,y in zip(aa,bb):
        if x!='-':pos+=1
        if x!=y:out.append({'Position in before (1-based)':str(pos) if x!='-' else f'after {pos}', 'Error':'Insertion' if x=='-' else 'Deletion' if y=='-' else 'Substitution','Before':x,'After':y})
    # Alignment positions describe net differences; repeated bases can admit ties.
    vu.alignment_view(r['Before'],r['After'])
    if out:st.dataframe(pd.DataFrame(out),hide_index=True,width='stretch')
    else:st.write('Sequences are identical.')
    with st.expander('Measurement details'):
        st.write('Top sequence = before; bottom = after. Positions use a minimum-edit alignment and describe net changes, not a unique mutation history. Up to four examples per strand/process are retained; event charts count the entire run. PCR examples show parent → child comparisons. Changed examples are prioritized.')

def stage_output(run,name):
    r=next(x for x in event_table(run) if x['Process']==name)
    a,b=st.columns(2);a.metric('New mutation events',r['sub']+r['ins']+r['del']);b.metric('Input nucleotides examined',f"{r['input_nt']:,}")
    chart_events([r],order=[name])
    st.write(f"S: {r['sub']:,} · I: {r['ins']:,} · D: {r['del']:,}")
    if name=='PCR':st.write('Denominator: bases copied into new children across all cycles.')
    elif name=='Sequencing':st.write('Denominator: template positions eligible for sequencing errors; R1/R2 counted separately.')
    else:st.write('Denominator: input bases in the selected error region, across all molecules.')

def art_controls(advanced,seed,key,compact=False):
    import art_sequencing as art
    paired=st.radio('Read mode',['Paired-end','Single-end'],horizontal=True,key=key+'_paired')=='Paired-end'
    length=st.number_input('Read length (nt)',25,150,125,key=key+'_length')
    if not compact:st.write('Quality / error model: ART empirical HS25 (HiSeq 2500)')
    options={}
    with st.expander('Advanced · ART configuration',expanded=advanced):
        executable=st.text_input('ART executable path (blank = PATH / ART_ILLUMINA)',key=key+'_executable')
        q1=st.number_input('R1 quality shift',-30,30,0,key=key+'_q1')
        q2=st.number_input('R2 quality shift',-30,30,0,key=key+'_q2',disabled=not paired)
        policy=st.selectbox('Short templates',['Shorten reads to template','Stop if template is short'],key=key+'_short')
        st.write('Shortening is an explicit variable-length approximation. It does not simulate adapter read-through; no bases are padded and no template is silently dropped.')
        timeout=st.number_input('ART total time budget (seconds)',1,3600,300,key=key+'_timeout')
        indel_policy=st.radio('ART sequencing indels',
                              ['Disabled — Windows-compatible default','Use ART native defaults','Set custom rates'],
                              horizontal=True,key=key+'_indel_policy')
        if indel_policy=='Use ART native defaults':
            options.update({'insertion_r1':None,'deletion_r1':None,'insertion_r2':None,'deletion_r2':None})
            st.warning('Some ART Windows 2.5.8 builds create invalid FASTQ when a sequencing indel changes read length. Use this only after checking your ART build.')
        elif indel_policy=='Set custom rates':
            for field,label,default in [('insertion_r1','R1 insertion',.009),('deletion_r1','R1 deletion',.011),
                                        ('insertion_r2','R2 insertion',.015),('deletion_r2','R2 deletion',.023)]:
                options[field]=st.number_input(label+' (%)',0.,10.,default,format='%.5f',key=key+'_'+field,disabled=field.endswith('r2') and not paired)/100
        else:
            st.write('ART still uses its empirical HS25 quality profile and creates substitutions from that profile. This setting disables only ART-generated sequencing insertions/deletions; DT4DDS molecular indels remain in the input template.')
    c=art.ArtSettings(executable=executable,seed=int(seed),paired=paired,read_length=int(length),q_shift_r1=int(q1),q_shift_r2=int(q2),
                      short_policy=policy,timeout_seconds=int(timeout),**options)
    found=art.detect_art_executable(executable)
    if found:
        if not compact:st.write('✓ ART executable found')
    else:st.info('ART executable not found. Select art_illumina.exe (Windows) or art_illumina (Linux/macOS). FASTQ upload remains available in Consensus.')
    return c


def results(rows,run,js,csv,zipped,current_public_settings):
    st.markdown('### Simulation results')
    if not run:
        return
    st.success('✓ Simulation complete')
    from read_pipeline import parse
    mates=[parse(run['r1'])]+([parse(run['r2'])] if run['r2'] else [])
    read_units=len(run.get('truth',[])) or len(mates[0])
    fastq_records=sum(map(len,mates))
    total_nt=sum(len(r['seq']) for m in mates for r in m)
    mean_q=sum(sum(r['q']) for m in mates for r in m)/total_nt if total_nt else 0
    mean_length=total_nt/fastq_records if fastq_records else 0

    designed_sids=[int(r['Strand ID']) for r in rows]
    unit_counts=pd.Series([int(t['sid']) for t in run.get('truth',[])],dtype='int64').value_counts()
    coverage=[int(unit_counts.get(sid,0)) for sid in designed_sids]
    represented=sum(x>0 for x in coverage)
    mean_coverage=sum(coverage)/len(coverage) if coverage else 0

    auto_evaluation_limit=10000
    if 'visual_evaluation' not in run and fastq_records<=auto_evaluation_limit:
        with st.spinner('Measuring FASTQ errors…'):
            run['visual_evaluation']=ev.sequencing(rows,run)
    report=run.get('visual_evaluation')
    totals=report['totals'] if report else None
    total_rate=(totals['errors']/totals['reference_nt']) if totals and totals['reference_nt'] else None

    a,b,c,d=st.columns(4)
    a.metric('Sampled strands',f'{represented} / {len(designed_sids)}')
    b.metric('Missing strands',len(designed_sids)-represented)
    c.metric('Mean coverage',f'{mean_coverage:.1f}×')
    d.metric('Error rate','Pending' if total_rate is None else f'{100*total_rate:.4f}%')
    st.write(
        f'{read_units:,} {"read pairs" if run["r2"] else "reads"} · '
        f'{fastq_records:,} FASTQ records · mean Q {mean_q:.2f} · mean length {mean_length:.2f} nt'
    )
    shortened=int(run.get('art',{}).get('shortened_read_units',0))
    if shortened:
        st.write(f'{shortened:,} read units used a shorter read length.')

    if report:
        records_by_sid={sid:[] for sid in designed_sids}
        for record in report['reads']:records_by_sid.setdefault(int(record['Strand ID']),[]).append(record)
        table=[]
        for sid,cov in zip(designed_sids,coverage):
            group=records_by_sid.get(sid,[])
            den=sum(r['reference_nt'] for r in group)
            sub=sum(r['S'] for r in group);ins=sum(r['I'] for r in group);dele=sum(r['D'] for r in group)
            errors=sum(r['errors'] for r in group)
            table.append({'Strand ID':sid,'Coverage':cov,'Substitution':sub,'Insertion':ins,
                          'Deletion':dele,'Total errors':errors,
                          'Error rate (%)':100*errors/den if den else None})
        st.dataframe(pd.DataFrame(table),hide_index=True,width='stretch',column_config={
            'Error rate (%)':st.column_config.NumberColumn(format='%.4f')})
    else:
        st.info(
            f'{fastq_records:,} FASTQ records exceed the automatic measurement limit.'
        )
        coverage_table=pd.DataFrame({'Strand ID':designed_sids,'Coverage':coverage})
        st.dataframe(coverage_table,hide_index=True,width='stretch')
        if st.button('Measure errors',key='art_measure_all'):
            with st.spinner('Measuring all records…'):
                run['visual_evaluation']=ev.sequencing(rows,run)
            st.rerun()

    with st.expander('View one read'):
        if report and report['reads']:
            inspect_sid=st.selectbox('Strand ID',designed_sids,key='result_inspect_sid')
            candidates=[r for r in report['reads'] if r['Strand ID']==inspect_sid]
            if candidates:
                selected=st.selectbox('FASTQ record',range(len(candidates)),
                    format_func=lambda i:candidates[i]['read_name'],key='result_inspect_record')
                record=candidates[selected]
                st.write(f'Mate {record["mate"]} · length {record["length"]} nt · mean Q {record["mean_Q"]:.2f} · S/I/D = {record["S"]}/{record["I"]}/{record["D"]}')
                vu.alignment_view(record['reference'],record['sequence'],record['comparison'].startswith('Observed'))
            else:st.warning('This strand received no read units.')
    public_manifest={'settings':run['settings'],'molecular_profile':run.get('preset'),
                     'art':{k:v for k,v in run['art'].items() if k not in ('templates','logs')}}
    files={'R1.fastq':run['r1'],'decoder_settings.json':js(asdict(current_public_settings(rows))),
           'evaluation/reference_strands.csv':csv(rows),'evaluation/expected_hash.txt':digest_bytes(run['stored_reference']).encode(),
           'channel_manifest.json':js(public_manifest),
           'evaluation/channel_config.json':js(public_manifest)}
    if run['r2']:files['R2.fastq']=run['r2']
    st.download_button('Download FASTQ package',zipped(files),'sequencing_output.zip',key='art_download')
    collect('ART run audit',js(run['art']),'art_run_audit.json')
    collect('Stage events',csv(run['stages']),'stage_events.csv')
    collect('Stage examples',csv(run['audit']),'stage_examples.csv')
    collect('PCR abundance',csv(run['abundance']),'pcr_abundance.csv')
    collect('Evaluation truth',js(run['truth']),'evaluation_truth.json')
    collect('R1 FASTQ',run['r1'],'R1.fastq')
    if run['r2']:collect('R2 FASTQ',run['r2'],'R2.fastq')


def digest_bytes(value):return hashlib.sha256(value).hexdigest()


def _legacy_render(rows,ui,digest,js,csv,zipped,current_public_settings):
    import channel_presets as presets
    import channel_runs
    with st.container(border=True):
        ui.step_header(5,'FASTQ Generation')
        if not rows:st.info('Prepare strands first.');return
        mode=st.radio('Channel model',['Illumina sequencing only','Cascaded DNA storage'],horizontal=True,key='channel_model')
        a,b=st.columns([3,1])
        with a:advanced=st.radio('Simulation controls',['Basic','Advanced'],horizontal=True,key='channel_controls')=='Advanced'
        with b:seed=st.number_input('Random seed',0,2**31-1,12345,key='channel_seed')
        design=digest(js(rows));signature=digest(js([design,mode,seed]))
        state=st.session_state.get('channel_steps')
        if not state or state.get('design')!=signature:
            state={'design':signature,'results':{},'signatures':{}};st.session_state.channel_steps=state
            st.session_state.pop('process_run',None);st.session_state.pop('sim_decode',None)
        box=st.empty();buttons=[];outputs=[];details=[]
        cfg=asdict(presets.defaults(int(seed)))
        if mode=='Illumina sequencing only':
            st.write('🧬 Clean oligos → 🔬 ART Illumina → 📄 FASTQ')
            st.write('Assumes error-free input molecules. Synthesis, storage and PCR mutations are off.')
            left,right=st.columns([1,1.35],gap='large')
            with left:
                artcfg=art_controls(advanced,seed,'illumina')
                coverage=st.number_input('Read pairs per oligo' if artcfg.paired else 'Reads per oligo',1,10000,350,key='illumina_coverage')
                st.write('Suggested levels: 5 · 10 · 20 · 50 · 100 · 350')
                clicked=st.button('Run Illumina Sequencing',type='primary',width='stretch',key='run_illumina')
            fingerprint=digest(js([asdict(artcfg),coverage]))
            if state['signatures'].get(0)!=fingerprint:
                state['results'].clear();st.session_state.pop('process_run',None);st.session_state.pop('sim_decode',None)
            if clicked:
                st.session_state.pop('process_run',None);st.session_state.pop('sim_decode',None)
                try:
                    with right:
                        bar=st.progress(0.,text='ART sequencing')
                        run=channel_runs.finish(rows,artcfg,reads_per_oligo=int(coverage),progress=bar.progress)
                    run.update(design_signature=design,stored_reference=st.session_state.get('stored_bytes',b''))
                    st.session_state.process_run=run;state['signatures'][0]=fingerprint;state['results'][0]=True
                    st.rerun()
                except Exception as exc:st.error(str(exc))
            with right:
                if st.session_state.get('process_run'):stage_output(st.session_state.process_run,'Sequencing')
                else:st.info('Choose coverage and read mode, then run ART.')
        else:
            st.selectbox('Paper preset',[presets.PRESET],key='paper_preset')
            with st.expander('Preset sources and assumptions',expanded=False):
                st.write('Gimpel et al., Nature Communications (2023). Simplified literature-derived scenario; not a full dt4dds reproduction or a lab calibration.')
                st.markdown('[Paper]('+presets.PAPER+') · [Reference implementation]('+presets.REPO+')')
                st.dataframe(pd.DataFrame(presets.provenance()),hide_index=True,width='stretch')
            cfg['scope']=st.radio('Molecular error region',['Full strand','BODY only'],horizontal=True,key='cascade_scope')
            if cfg['scope']=='BODY only':st.info('BODY-only is a controlled assumption: SI and primers are protected from molecular errors. ART still reads the full template and may introduce errors anywhere.')
            if cfg['scope']=='Full strand' and (st.session_state.get('hedges_meta') or {}).get('sid_mode')==hb.SID_MODE_NO_HAMMING:
                st.warning('The unprotected SI can receive molecular and sequencing errors.')
            # The Basic view intentionally shows only two actions: make the
            # molecular pool, then sequence it.  The exact same Settings and
            # checkpoints are used by the Advanced stage-by-stage view below.
            if not advanced:
                st.write('🧬 Synthesis  →  🧊 Storage  →  🧪 PCR  →  🔬 ART sequencing')
                st.caption('Literature-derived preset. Detailed rates and per-process inspection are in Advanced.')
                cfg['founders']=st.select_slider('Starting molecules per strand',options=[1,2,3],value=2,key='basic_founders')
                molecular_key=digest(js(asdict(pm.Settings(**cfg))))
                if state['signatures'].get('molecular')!=molecular_key:
                    state['results'].clear();state['signatures'].clear()
                    st.session_state.pop('process_run',None);st.session_state.pop('sim_decode',None)
                left,right=st.columns([1,1.35],gap='large')
                with left:
                    st.metric('PCR cycles',cfg['cycles'])
                    st.metric('Storage exposure','Off' if not cfg['storage'] else f"{cfg['storage_halflives']} half-lives")
                    make_pool=st.button('Run molecular processes',type='primary',width='stretch',key='run_molecular_processes')
                with right:
                    if make_pool:
                        try:
                            bar=st.progress(0.,text='Building molecular pool')
                            checkpoint=pm.simulate(rows,pm.Settings(**cfg),stop_after='PCR',stage_progress=lambda n,f:bar.progress(f,text=n))
                            state['results'][2]=checkpoint;state['signatures']['molecular']=molecular_key
                            st.rerun()
                        except Exception as exc:st.error(str(exc))
                    checkpoint=state['results'].get(2)
                    if checkpoint:
                        st.success('✓ Molecular pool ready')
                        compact=pd.DataFrame(event_table(checkpoint)).set_index('Process')[['sub','ins','del']]
                        st.bar_chart(compact)
                    else:st.info('Run the molecular processes once to create the PCR pool.')
                st.markdown('### 5.4 ART sequencing')
                seq_left,seq_right=st.columns([1,1.35],gap='large')
                with seq_left:
                    artcfg=art_controls(False,seed,'cascade_art_basic')
                    sampling_mode=st.radio('Sequencing sampling',['Equal per strand','Pool-proportional'],horizontal=True,key='basic_sampling')
                    coverage=st.number_input('Reads per strand' if not artcfg.paired else 'Read pairs per strand',1,10000,350,key='basic_coverage')
                    total=int(coverage)*len(rows)
                    st.write(f'App will generate {total:,} '+('read pairs' if artcfg.paired else 'reads')+f' for {len(rows)} designed strands.')
                    if sampling_mode=='Equal per strand':st.caption('Every surviving SID receives the same number of read units; PCR variants within that SID remain abundance-weighted.')
                    else:st.caption('Read units are drawn from the whole PCR pool; high-abundance SIDs can receive more reads and low-abundance SIDs can be absent.')
                    run_seq=st.button('Run sequencing',type='primary',width='stretch',key='run_cascade_basic_sequencing',disabled=not bool(checkpoint))
                with seq_right:
                    if run_seq and checkpoint:
                        try:
                            bar=st.progress(0.,text='Grouping variable-length templates and running ART')
                            run=channel_runs.finish(rows,artcfg,checkpoint=checkpoint,total_units=total,progress=bar.progress,sampling_mode=sampling_mode)
                            run.update(design_signature=design,stored_reference=st.session_state.get('stored_bytes',b''),preset=presets.manifest(pm.Settings(**cfg)))
                            state['results'][3]=run;state['signatures']['sequencing']=digest(js([asdict(artcfg),total,molecular_key]));st.session_state.process_run=run
                            st.rerun()
                        except Exception as exc:st.error(str(exc))
                    if st.session_state.get('process_run'):
                        stage_output(st.session_state.process_run,'Sequencing')
                    elif checkpoint:st.info('The pool is ready. Choose read mode and run ART sequencing.')
                states=['✓ Done' if checkpoint else 'Ready','✓ Done' if checkpoint else 'Ready','✓ Done' if checkpoint else 'Ready','✓ Done' if st.session_state.get('process_run') else ('Ready' if checkpoint else 'Waiting')]
                pipeline(box,states)
                results(rows,st.session_state.get('process_run'),js,csv,zipped,current_public_settings)
                return
            def rate(label,k):cfg[k]=st.number_input(label+' (%)',0.,100.,float(cfg[k]*100),format='%.7f',key='paper_'+k)/100
            for i,name in enumerate(STAGES):
                st.markdown(f'### 5.{i+1} '+('Sampling → Illumina sequencing' if i==3 else name))
                left,right=st.columns([1,1.35],gap='large')
                with left:
                    if i==0:
                        cfg['founders']=st.number_input('Mean starting molecules per strand',1,100,20,key='paper_founders')
                        cfg['abundance_sigma']=st.number_input('Abundance variation (lognormal sigma)',0.,3.,.30,format='%.3f',key='paper_abundance_sigma')
                        st.write('Fixed molecule budget; one starting molecule is reserved per SID. This is a no-synthesis-dropout approximation.')
                        rate('Synthesis deletion','syn_del')
                        with st.expander('Advanced · synthesis',expanded=advanced):
                            cfg['synthesis']=st.checkbox('Enable synthesis defects',True,key='paper_synthesis')
                            rate('Synthesis substitution','syn_sub');rate('Synthesis insertion','syn_ins')
                            st.write('Independent-site approximation; deletion-run and sequence-context biases from dt4dds are not reproduced.')
                    elif i==1:
                        cfg['storage']=st.checkbox('Apply storage exposure',False,key='paper_storage')
                        if cfg['storage']:
                            cfg['storage_halflives']=st.number_input('Storage exposure (half-lives)',0.,10.,1.,step=.5,key='paper_halflives')
                            loss=st.checkbox('Include molecule loss',False,key='paper_storage_loss')
                            cfg['storage_survival']=2**(-cfg['storage_halflives']) if loss else 1.
                            st.write(f"Literature-derived substitution probability: {100*cfg['storage_rate_per_halflife']*cfg['storage_halflives']:.5f}%. Survival: {100*cfg['storage_survival']:.2f}%.")
                            if not loss:st.write('Sequence-damage study: all molecules retained. No claim of physical survival at this exposure.')
                            with st.expander('Advanced · storage',expanded=advanced):
                                rate('Substitution per half-life','storage_rate_per_halflife')
                                rate('Additional C → T per C base over exposure','storage_ct')
                                rate('Additional substitution per base over exposure','storage_sub')
                                st.write('Uniform substitution approximation to the source model. No calendar-time or temperature calibration; loss removes whole molecules, not explicit fragments.')
                        else:st.write('No storage exposure assumed; molecules pass through unchanged.')
                    elif i==2:
                        cfg['cycles']=st.number_input('PCR cycles',0,20,10,key='paper_cycles')
                        cfg['efficiency']=st.number_input('Copy probability per molecule per cycle (%)',0.,100.,80.,key='paper_efficiency')/100
                        cfg['efficiency_relative_sd']=st.number_input('Relative amplification-factor SD',0.,.5,.0051,format='%.4f',key='paper_efficiency_sd')
                        st.write('Q5-derived mutation rate; parent errors are inherited by descendants.')
                        with st.expander('Advanced · PCR',expanded=advanced):
                            cfg['pcr']=st.checkbox('Enable PCR amplification',True,key='paper_pcr')
                            rate('Substitution per newly copied base','pcr_sub');rate('Insertion per copied position','pcr_ins');rate('Deletion per copied base','pcr_del')
                            cfg['max_variants']=st.number_input('Variant / work budget',1000,1000000,100000,step=1000,key='paper_budget')
                    else:
                        artcfg=art_controls(advanced,seed,'cascade_art')
                        total=st.number_input('Total read pairs' if artcfg.paired else 'Total reads',1,100000,min(100000,350*len(rows)),key='paper_total')
                        st.write(f'Target mean: {total/len(rows):.2f} per designed oligo. Actual coverage follows the pool abundance.')
                    buttons.append(st.button('Run '+name,type='primary',width='stretch',key='run_channel_'+name))
                with right:outputs.append(st.empty())
                details.append(st.empty())
            fields=[['seed','scope','founders','abundance_sigma','synthesis','syn_sub','syn_ins','syn_del'],
                    ['storage','storage_sub','storage_ct','storage_survival','storage_halflives','storage_rate_per_halflife'],
                    ['pcr','cycles','efficiency','efficiency_relative_sd','pcr_sub','pcr_ins','pcr_del','max_variants']]
            fingerprints=[digest(js({k:cfg[k] for k in f})) for f in fields]+[digest(js([asdict(artcfg),total]))]
            for i in range(4):
                if i in state['results'] and state['signatures'].get(i)!=fingerprints[i]:
                    for j in range(i,4):state['results'].pop(j,None);state['signatures'].pop(j,None)
                    st.session_state.pop('process_run',None);st.session_state.pop('sim_decode',None);break
            states=['✓ Done' if i in state['results'] else 'Ready' if i==0 or i-1 in state['results'] else 'Waiting' for i in range(4)]
            pipeline(box,states)
            cfgobj=pm.Settings(**cfg)
            source=presets.manifest(cfgobj)
            st.write('Parameter status: '+('Custom — based on paper preset' if source['overrides'] else 'Literature-derived preset (simplified)'))
            for i,clicked in enumerate(buttons):
                if not clicked:continue
                if i and i-1 not in state['results']:st.warning('Run '+STAGES[i-1]+' first.');continue
                for j in range(i,4):state['results'].pop(j,None);state['signatures'].pop(j,None)
                st.session_state.pop('process_run',None);st.session_state.pop('sim_decode',None)
                states[i]='⏳ Running';pipeline(box,states)
                try:
                    with outputs[i].container():
                        bar=st.progress(0.,text=STAGES[i])
                        if i<3:
                            run=pm.simulate(rows,cfgobj,stop_after=STAGES[i],checkpoint=state['results'].get(i-1),stage_progress=lambda n,f:bar.progress(f,text=n))
                        else:run=channel_runs.finish(rows,artcfg,checkpoint=state['results'][2],total_units=int(total),progress=bar.progress)
                    state['results'][i]=run;state['signatures'][i]=fingerprints[i]
                    if i==3:
                        run.update(design_signature=design,stored_reference=st.session_state.get('stored_bytes',b''),preset=source)
                        st.session_state.process_run=run
                    st.rerun()
                except Exception as exc:
                    states[i]='Error';pipeline(box,states);st.error(str(exc))
            for i,name in enumerate(STAGES):
                run=state['results'].get(i)
                with outputs[i].container():
                    if run:
                        stage_output(run,name)
                        if i==1:
                            lost=sum(r.get('lost_molecules',0) for r in run['stages'] if r['Stage']=='Storage')
                            st.metric('Molecules lost',lost)
                    else:st.info('Run '+name+' to view results.')
                if run:
                    with details[i].container():before_after(run,name)
        results(rows,st.session_state.get('process_run'),js,csv,zipped,current_public_settings)


def render(rows,ui,digest,js,csv,zipped,current_public_settings):
    """Single defensible workflow: official DT4DDS molecular channel + real ART."""
    import channel_runs
    import dt4dds_backend as molecular

    with st.container(border=True):
        ui.step_header(5,'FASTQ Generation')
        if not rows:
            st.info('Prepare strands first.')
            return

        cfg=molecular.MolecularSettings()
        design=digest(js(rows))
        st.write('Model: published experimental DNA-storage profile')

        control_left,control_right=st.columns(2,gap='large')
        with control_left:
            coverage=st.number_input('Target coverage per strand',1,10000,350,key='experimental_coverage')
            sampling_label=st.radio('Sampling method',['Pool distribution','Equal coverage'],horizontal=True,
                key='experimental_sampling',help='Pool distribution follows molecular abundance. Equal coverage assigns the same target read units to each surviving strand.')
            with st.expander('Advanced · reproducibility',expanded=False):
                seed=st.number_input('Random seed',0,2**31-1,cfg.seed,key='experimental_seed')
            cfg=molecular.MolecularSettings(seed=int(seed))
        with control_right:
            artcfg=art_controls(False,cfg.seed,'experimental_art',compact=True)
            artcfg=replace(artcfg,seed=int(seed))
        total=int(coverage)*len(rows)
        sampling_mode='Pool-proportional' if sampling_label=='Pool distribution' else 'Equal per strand'
        ready,dependency_message=molecular.dependency_status()
        art_ready=__import__('art_sequencing').detect_art_executable(artcfg.executable)
        if not ready:st.error(dependency_message)
        run_clicked=st.button('Run simulation',type='primary',width='stretch',key='run_experimental_channel',disabled=not (ready and art_ready))

        signature=digest(js([design,asdict(cfg),asdict(artcfg),total,sampling_mode]))
        if st.session_state.get('channel_signature')!=signature:
            st.session_state.channel_signature=signature
            st.session_state.pop('process_run',None)
            st.session_state.pop('sim_decode',None)

        if run_clicked:
            st.session_state.pop('process_run',None)
            st.session_state.pop('sim_decode',None)
            try:
                bar=st.progress(0.,text='Running simulation')
                checkpoint=molecular.run(rows,cfg,progress=lambda _name,f:bar.progress(f,text='Running simulation'))
                bar.progress(.82,text='Generating FASTQ')
                run=channel_runs.finish(rows,artcfg,checkpoint=checkpoint,total_units=total,
                    progress=lambda f:bar.progress(.82+.18*f,text='Generating FASTQ'),
                    sampling_mode=sampling_mode)
                run.update(design_signature=design,stored_reference=st.session_state.get('stored_bytes',b''))
                st.session_state.process_run=run
                st.rerun()
            except Exception as exc:
                st.error(str(exc))

        results(rows,st.session_state.get('process_run'),js,csv,zipped,current_public_settings)
