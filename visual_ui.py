from exports import collect
"""Visual inspection and post-hoc evaluation; no reference enters decoding."""
import html
import json
import pandas as pd
import streamlit as st
import evaluation as ev


def flow(settings=None,done=False):
    settings=settings or {}
    cards=[('🧬','1 · Synthesis','synthesis','Create molecules'),('🧊','2 · Storage','storage','Modify stored molecules'),
           ('🧪','3 · PCR','pcr','Copy parent sequences'),('🔬','4 · Sequencing','sequencing','Generate FASTQ reads')]
    st.markdown('''<style>.dna-flow{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:12px 0 20px}
    .dna-stage{border:1px solid #b9cbdc;border-top:4px solid #328aab;border-radius:12px;padding:16px;background:var(--secondary-background-color)}
    .dna-stage .icon{font-size:28px}.dna-stage p{margin:8px 0}.dna-stage small{opacity:.8}</style>''',unsafe_allow_html=True)
    blocks=[]
    for icon,name,key,desc in cards:
        on=settings.get(key,key!='storage')
        status=('✓ Completed' if done else '● Enabled') if on else ('✓ Passed through · errors off' if done else '○ Errors off / PCR bypassed')
        blocks.append(f'<div class="dna-stage"><div class="icon">{icon}</div><b>{name}</b><p>{status}</p><small>{desc}</small></div>')
    st.markdown('<div class="dna-flow">'+''.join(blocks)+'</div>',unsafe_allow_html=True)


def alignment_view(reference,observed,partial=False):
    a,b=ev.alignment(reference,observed,partial)
    colors={'S':'#ffccd2','I':'#bce8fa','D':'#ffe4a3','N':'#e1ccff'}
    top=[];bottom=[]
    st.write('🟥 Substitution · 🟦 Insertion · 🟨 Deletion · 🟪 Unresolved N')
    for x,y in zip(a,b):
        kind='N' if y=='N' else 'I' if x=='-' else 'D' if y=='-' else 'S' if x!=y else ''
        style=f'background:{colors[kind]};color:#17212b;' if kind else ''
        top.append(f'<span style="{style}">{html.escape(x)}</span>')
        bottom.append(f'<span style="{style}">{html.escape(y)}</span>')
    st.markdown('<div style="overflow-x:auto;font-family:monospace;white-space:pre;line-height:1.9;padding:12px;border:1px solid #b9cbdc;border-radius:8px">Reference  '+''.join(top)+'<br>Observed   '+''.join(bottom)+'</div>',unsafe_allow_html=True)


def sequence_results(rows,run):
    st.markdown('#### 🔬 Sequencing output · inspect the reads')
    if 'visual_evaluation' not in run:
        with st.spinner('Measuring read errors against references (evaluation only)…'):
            run['visual_evaluation']=ev.sequencing(rows,run)
    report=run['visual_evaluation'];totals=report['totals']
    st.write('Final read errors, including inherited errors. R1 and R2 are counted separately.')
    cols=st.columns(4)
    for col,key,label in zip(cols,('S','I','D','errors'),('Substitutions (nt)','Insertions (nt)','Deletions (nt)','Total read edits')):col.metric(label,totals[key])
    strict=[r['Strict similarity (%)'] for r in report['reads']]
    st.metric('Before consensus · total edit rate',f"{100*totals['errors']/totals['reference_nt']:.4f}%" if totals['reference_nt'] else 'N/A')
    with st.expander('Advanced · strict positional similarity'):
        st.metric('Strict read similarity',f'{sum(strict)/len(strict):.2f}%' if strict else 'N/A')
        st.write('Different length = 0%; otherwise positional matches without alignment. This is not an edit-error rate.')
    a,b=st.columns(2)
    with a:
        st.markdown('**Coverage distribution**')
        coverage_counts=pd.Series([r['Read units'] for r in report['coverage']]).value_counts().sort_index()
        st.bar_chart(coverage_counts.rename_axis('Read units per strand').rename('Number of strands'))
    with b:
        st.markdown('**Final read errors by strand**')
        st.bar_chart(pd.DataFrame(report['reads']).groupby('Strand ID')[['S','I','D']].sum())
    a,b=st.columns(2)
    with a:
        st.markdown('**Coverage by strand**')
        frame=pd.DataFrame(report['coverage'])
        st.bar_chart(frame.set_index('Strand ID')[['Read units']])
    with b:
        st.dataframe(frame,hide_index=True,width='stretch')
    sid=st.selectbox('Inspect sequencing strand',list(frame['Strand ID']),key='seq_inspect_sid')
    reads=[r for r in report['reads'] if r['Strand ID']==sid]
    if not reads:
        st.warning('This strand received zero sampled reads.');return
    st.dataframe(pd.DataFrame([{k:v for k,v in r.items() if k!='reference'} for r in reads]),hide_index=True,width='stretch')
    selected=st.selectbox('Read alignment',range(len(reads)),format_func=lambda i:reads[i]['read_name'],key='seq_inspect_read')
    row=reads[selected];alignment_view(row['reference'],row['sequence'],row['comparison'].startswith('Observed'))
    collect('Download all read measurements',pd.DataFrame(report['reads']).drop(columns='reference').to_csv(index=False).encode(),'read_errors.csv')
    fastq_preview(run['r1'],run['r2'])


def consensus_results(result,prefix,rows=None):
    st.markdown('#### 🧩 Consensus · one sequence per strand')
    frame=pd.DataFrame(result.get('consensus_rows',[]))
    if not frame.empty:
        frame['Sequence']=[result['all_consensus'].get(int(s),'') for s in frame['Strand ID']]
        st.dataframe(frame,hide_index=True,width='stretch')
    if rows is None:
        st.info('Experimental FASTQ: true S/I/D counts and percentage error reduction require a known reference. Support, N counts and decoder status are available without one.')
        return
    report=ev.consensus_report(rows,result)
    status_counts=pd.Series([r['Outcome'] for r in report['strands']]).value_counts()
    badges=st.columns(4)
    for col,status in zip(badges,('EXACT','ERRORS REMAIN','UNCERTAIN','NO CONSENSUS')):
        col.metric(status.title(),int(status_counts.get(status,0)))
    a,b,c=st.columns(3)
    fmt=lambda x:'N/A' if x is None else f'{100*x:.3f}%'
    a.metric('Before consensus · BODY error rate',fmt(report['before_rate']))
    b.metric('After consensus · BODY error rate',fmt(report['after_rate']))
    c.metric('Error-rate reduction','N/A' if report['reduction_percent'] is None else f'{report["reduction_percent"]:.2f}%')
    st.write(f'Compared: {report["matched_strands"]}/{report["total_strands"]} strands. Missing consensus is a failure, not zero errors.')
    exact=int(status_counts.get('EXACT',0))
    st.metric('Exact consensus strands',f'{100*exact/report["total_strands"]:.1f}%' if report['total_strands'] else 'N/A')
    if report['matched_strands']:
        st.bar_chart(pd.DataFrame({'Error rate (%)':[100*report['before_rate'],100*report['after_rate']]},index=['Before consensus','After consensus']))
    st.dataframe(pd.DataFrame(report['strands']),hide_index=True,width='stretch')
    with st.expander('Error counts before and after consensus · S / I / D / N',expanded=True):
        counts=[]
        for k,label in [('S','Substitution'),('I','Insertion'),('D','Deletion'),('N','Unresolved N')]:
            counts.append({'Type':label,
                           'Before: all assigned reads':sum(r['Before '+k] for r in report['strands']),
                           'After: nonempty consensus':sum(r['After '+k] or 0 for r in report['strands'])})
        st.dataframe(pd.DataFrame(counts),hide_index=True,width='stretch')
        st.write('Many reads → one consensus: compare normalized rates above, not the raw count ratio.')
    report['assignment_audit_note']='Sequencing truth and cluster assignment are evaluated separately; references never select a consensus.'
    collect('Download consensus evaluation',json.dumps(report,indent=2).encode(),'consensus_evaluation.json',key=prefix+'_eval_download')
    document='<!doctype html><meta charset="utf-8"><title>Consensus evaluation</title><style>body{font:15px system-ui;margin:35px}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #bbb;padding:6px}</style><h1>Consensus evaluation</h1>'
    document+='<p>'+html.escape(report['definition'])+'</p>'
    document+='<pre>'+html.escape(json.dumps({k:v for k,v in report.items() if k!='strands'},indent=2))+'</pre>'
    document+=pd.DataFrame(report['strands']).drop(columns='Consensus').to_html(index=False)
    document+='<h2>Public decoder settings</h2><pre>'+html.escape(json.dumps(result['settings'],indent=2))+'</pre>'
    collect('Presentation report',document.encode(),'consensus_report.html')
    with st.expander('Advanced · metric definition and alignment',expanded=False):
        st.write(report['definition'])
        refs={int(r['Strand ID']):str(r['Payload'])+str(r.get('Filler','')) for r in rows}
        sid=st.selectbox('Inspect consensus strand',list(refs),key=prefix+'_consensus_sid')
        seq=result['all_consensus'].get(sid,'')
        if seq:alignment_view(refs[sid],seq)
        else:st.warning('No consensus for this strand. It is not counted as a zero-error success.')


def fastq_preview(r1,r2=b''):
    with st.expander('FASTQ preview · quality scores'):
        mate=st.radio('Preview mate',['R1','R2'] if r2 else ['R1'],horizontal=True)
        records=ev.rp.parse(r2 if mate=='R2' else r1)
        if not records:return
        record=records[0]
        st.write('Four lines: read name / DNA sequence / + / quality characters.')
        st.code('@'+record['name']+'\n'+record['seq']+'\n+\n'+''.join(chr(q+33) for q in record['q']),language=None)
        st.write('Numeric Q = character code − 33. One quality value per nucleotide.')
        st.dataframe(pd.DataFrame({'Position':range(1,len(record['seq'])+1),'Base':list(record['seq']),'Q':record['q']}),hide_index=True,width='stretch')


def length_summary(lengths,label):
    if not lengths:return
    a,b,c,d=st.columns(4)
    a.metric('Sequences',len(lengths));b.metric('Min length',min(lengths));c.metric('Mean length',f'{sum(lengths)/len(lengths):.1f}');d.metric('Max length',max(lengths))
    with st.expander(label):
        st.bar_chart(pd.Series(lengths).value_counts().sort_index().rename_axis('Length (nt)').rename('Count'))


def length_comparison(result):
    table=[]
    for sid,seq in result.get('all_consensus',{}).items():
        lengths=[len(r['sequence']) for r in result['observed_bodies'] if r['sid']==sid]
        table.append({'Strand ID':sid,'Reads before':len(lengths),
                      'Before min nt':min(lengths) if lengths else None,
                      'Before mean nt':sum(lengths)/len(lengths) if lengths else None,
                      'Before max nt':max(lengths) if lengths else None,
                      'After consensus nt':len(seq) if seq else None})
    st.markdown('#### Length before / after consensus · BODY only')
    st.dataframe(pd.DataFrame(table),hide_index=True,width='stretch')
    if table:
        st.line_chart(pd.DataFrame(table).set_index('Strand ID')[['Before mean nt','After consensus nt']])
    collect('Length comparison',pd.DataFrame(table).to_csv(index=False).encode(),'length_comparison.csv')


def stage_results(run):
    st.markdown('#### 🧬 Synthesis / 🧊 Storage · before and after')
    st.write('Measured region: '+run['settings']['scope'])
    data=run.get('stage_quality',[])
    for stage in ('Synthesis','Storage'):
        group=[r for r in data if r['Stage']==stage]
        if not group:continue
        st.markdown('**'+stage+(' · completed' if group[0]['Enabled'] else ' · disabled')+'**')
        rows=[]
        for side in ('Before','After'):
            count=sum(r[side+' molecules'] for r in group);den=sum(r[side+' reference_nt'] for r in group)
            errors=sum(r.get(side+' errors',0) for r in group)
            rows.append({'State':side,'S':sum(r.get(side+' S',0) for r in group),
                         'I':sum(r.get(side+' I',0) for r in group),'D':sum(r.get(side+' D',0) for r in group),
                         'Error rate (%)':100*errors/den if den else 0,
                         'Exact molecules (%)':sum(r[side+' exact_molecules_percent']*r[side+' molecules'] for r in group)/count})
        st.dataframe(pd.DataFrame(rows),hide_index=True,width='stretch')
    st.write('Before/after compare molecules with the designed sequence. New-event counts below identify where errors were introduced.')
    with st.expander('Per-strand stage measurements'):
        st.dataframe(pd.DataFrame(data),hide_index=True,width='stretch')


def presentation_report(rows,run,result):
    """Standalone, printable evidence with saved settings, not current widget edits."""
    report=ev.consensus_report(rows,result)
    table=pd.DataFrame(report['strands']).drop(columns=['Consensus'])
    payloads=ev.hedges_report(rows,result) if result.get('strand_results') else []
    summary={k:report[k] for k in ('matched_strands','total_strands','before_rate','after_rate','reduction_percent')}
    exact=result.get('data')==run.get('stored_reference') if result.get('data') is not None else None
    raw=run['visual_evaluation']['totals']
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>DNA storage experiment</title>
    <style>body{{font:15px system-ui;max-width:1200px;margin:40px auto;padding:0 20px;color:#152d3c}}h1,h2{{color:#166080}}table{{border-collapse:collapse;font-size:12px;display:block;overflow:auto}}th,td{{border:1px solid #ccd9df;padding:7px;text-align:left}}th{{background:#e9f3f7}}pre{{background:#f1f6f8;padding:15px;white-space:pre-wrap}}@media print{{body{{margin:0}}table{{font-size:9px}}}}</style>
    <h1>DNA storage · measured simulation run</h1>
    <p>Synthesis → storage → PCR → sequencing → consensus → HEDGES/file decoding</p>
    <p>Synthetic controlled coverage, no strand dropout or outer ECC. This single run does not establish wet-lab accuracy or a statistical recovery guarantee.</p>
    <h2>1. Final FASTQ errors</h2><pre>{html.escape(json.dumps(raw,indent=2))}</pre>
    <p>Minimum-edit alignment counts across all FASTQ records; R1/R2 counted separately. Short reads use observed-span alignment.</p>
    <h2>2. Consensus reduction</h2><pre>{html.escape(json.dumps(summary,indent=2))}</pre>
    <p>Rates above are fractions. {html.escape(report['definition'])}</p>{table.to_html(index=False)}
    <h2>3. File and payload checks</h2><p>Stored bytes exactly match: {exact if exact is not None else 'Not reconstructed'}. Decoder status: {html.escape(result['status'])}.</p>
    {pd.DataFrame(payloads).to_html(index=False) if payloads else '<p>HEDGES per-strand payload results unavailable.</p>'}
    <h2>4. New mutation events by process</h2>{pd.DataFrame(run['stages']).fillna(0).to_html(index=False)}
    <p>Event counts are not final read errors: descendants inherit upstream mutations; different event histories can produce the same final sequence.</p>
    <h2>5. Reproducible settings</h2><h3>Simulation</h3><pre>{html.escape(json.dumps(run['settings'],indent=2))}</pre>
    <h3>Public decoder settings</h3><pre>{html.escape(json.dumps(result['settings'],indent=2))}</pre>
    <p>Ground truth was used only for post-hoc evaluation, never for read assignment, consensus selection or file reconstruction.</p></html>'''
    return body.encode()


def strict_consensus_summary(result,rows):
    if rows is None:
        st.info('Upload reference strands to measure similarity. Consensus support is available without a reference.')
        return
    table=[]
    for row in rows:
        sid=int(row['Strand ID']);ref=str(row['Payload'])+str(row.get('Filler',''))
        reads=[r['sequence'] for r in result.get('observed_bodies',[]) if r['sid']==sid]
        consensus=result.get('all_consensus',{}).get(sid,'')
        table.append({'Strand ID':sid,'Before (%)':sum(ev.strict_similarity(ref,x) for x in reads)/len(reads) if reads else 0.,
                      'After (%)':ev.strict_similarity(ref,consensus) if consensus else 0.,'Reads':len(reads),'Consensus available':bool(consensus)})
    if not table:return
    before=sum(r['Before (%)'] for r in table)/len(table);after=sum(r['After (%)'] for r in table)/len(table)
    a,b,c=st.columns(3)
    a.metric('Before · strict BODY similarity',f'{before:.2f}%')
    b.metric('After · strict BODY similarity',f'{after:.2f}%')
    c.metric('Improvement',f'{after-before:+.2f} pp')
    import altair as alt
    frame=pd.DataFrame({'Step':['Before consensus','After consensus'],'Similarity (%)':[before,after]})
    st.altair_chart(alt.Chart(frame).mark_bar().encode(x=alt.X('Step:N',sort=['Before consensus','After consensus']),y=alt.Y('Similarity (%):Q',scale=alt.Scale(domain=[0,100])),tooltip=['Step','Similarity (%)']),use_container_width=True)
    st.write('Same BODY region, averaged equally across designed strands. Different length = 0%; missing reads/consensus = 0%. Improvement is in percentage points (pp). Step 5.5 measures full read spans, so its score has a different denominator.')
    st.dataframe(pd.DataFrame(table),hide_index=True,width='stretch')


def consensus_summary(result,rows):
    if rows is None:
        st.info('A reference is needed to measure true error reduction. Decoder processing never uses this reference.');return
    report=ev.consensus_report(rows,result)
    a,b,c=st.columns(3)
    fmt=lambda x:'N/A' if x is None else f'{100*x:.4f}%'
    a.metric('Before · BODY edit rate',fmt(report['before_rate']))
    b.metric('After · BODY edit rate',fmt(report['after_rate']))
    c.metric('Error-rate reduction','N/A' if report['reduction_percent'] is None else f"{report['reduction_percent']:.2f}%")
    exact=sum(x['Outcome']=='EXACT' for x in report['strands'])
    missing=sum(x['Outcome']=='NO CONSENSUS' for x in report['strands'])
    a,b=st.columns(2);a.metric('Exact consensus / designed strands',f"{exact} / {report['total_strands']}");b.metric('Missing consensus',missing)
    if report['before_rate'] is not None:
        st.bar_chart(pd.DataFrame({'BODY edit rate (%)':[100*report['before_rate'],100*report['after_rate']]},index=['Before','After']))
    st.write(f"Rates compare the same {report['matched_strands']} nonempty SID groups; missing groups are reported separately.")
    with st.expander('Advanced · strict positional similarity'):
        strict_consensus_summary(result,rows)
