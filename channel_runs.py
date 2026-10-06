"""Join process checkpoints to ART, retaining the existing app result contract."""
from collections import Counter
from dataclasses import asdict
import copy
import numpy as np
import art_sequencing as art


def finish(rows,cfg,reads_per_oligo=350,checkpoint=None,total_units=None,progress=None,sampling_mode='Pool-proportional'):
    if checkpoint is None:
        templates=[{'sid':int(r['Strand ID']),'template':r['Full strand'],'count':int(reads_per_oligo)} for r in rows]
        run={'settings':{'scope':'Full strand'},'stages':[],'stage_quality':[],'abundance':[
            {'Strand ID':int(r['Strand ID']),'Initial molecules':None,'Final molecules':None} for r in rows],
            'audit':[],'trace':[]}
        mode='Illumina sequencing only'
    else:
        if checkpoint.get('engine')=='DT4DDS':
            if checkpoint.get('completed_stage')!='Sampling':raise ValueError('Complete DT4DDS physical sampling before sequencing')
            source_templates=checkpoint.get('_templates',[])
        else:
            if checkpoint.get('completed_stage')!='PCR':raise ValueError('Complete PCR before sequencing')
            source_templates=art.pool_templates(rows,checkpoint)
        # Separate deterministic random stream: changing ART quality does not
        # change the sampled molecular population.
        seed=int(np.random.SeedSequence([cfg.seed,71]).generate_state(1)[0])
        templates=art.allocate(source_templates,int(total_units),seed,sampling_mode)
        run={k:copy.deepcopy(v) for k,v in checkpoint.items() if not k.startswith('_')}
        mode=checkpoint.get('settings',{}).get('mode','Cascaded DNA storage')
    output=art.run_art(templates,cfg,progress)
    run.update(output)
    run['stages'].append(output['sequencing_stage'])
    for unused in ('seq_sub','seq_ins','seq_del','quality_mode','q_profile','q_fixed','q_start','q_end','q_std'):
        run['settings'].pop(unused,None)
    run['settings'].update({'sequencing_engine':'ART','mode':mode,'paired':cfg.paired,'read_length':cfg.read_length,
                            'seed':cfg.seed,'sampling':sampling_mode if checkpoint else 'Fixed coverage',
                            'coverage':len(output['truth'])/len(rows),'art':asdict(cfg)})
    counts=Counter(t['sid'] for t in output['truth'])
    for row in run['abundance']:row['Read units']=counts[row['Strand ID']]
    audit_counts=Counter()
    for t in output['truth']:
        for mate,(span,errors) in enumerate(zip(t['sequencing_reference'],t['sequencing_edits'])):
            if audit_counts[t['sid']]>=4:break
            audit_counts[t['sid']]+=1
            # Bounded examples only; aggregate events come from every SAM record.
            run['audit'].append({'Stage':'Sequencing','Strand ID':t['sid'],'Founder':'sampled variant',
                'Cycle':None,'Before':span,'After':'','Before nt':len(span),'After nt':0,
                'New S':errors['S'],'New I':errors['I'],'New D':errors['D'],'unit':t['unit'],'mate':mate})
    from read_pipeline import parse
    mates=[parse(run['r1'])]+([parse(run['r2'])] if cfg.paired else [])
    for r in run['audit']:
        if r['Stage']=='Sequencing':
            r['After']=mates[r['mate']][r['unit']]['seq'];r['After nt']=len(r['After'])
    run['completed_stage']='Sequencing'
    return run
