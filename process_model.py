"""Sequence-error cascade with optional strand-specific PCR efficiency.

Every oligo has founders; PCR retains parents. Fixed coverage and pool sampling
are explicit modes; the latter can produce zero-read strands. No physical
strand dropout or outer ECC is introduced. Rates are not vendor-calibrated.
"""
from dataclasses import dataclass, asdict
from collections import Counter
import math
import copy
import numpy as np


@dataclass(frozen=True)
class Settings:
    seed: int = 7
    founders: int = 5
    coverage: int = 10
    scope: str = "BODY only"
    synthesis: bool = True
    syn_sub: float = 0.0001
    syn_ins: float = 0.0001
    syn_del: float = 0.001
    storage: bool = False
    storage_sub: float = 0.0001
    storage_ct: float = 0.001
    pcr: bool = True
    cycles: int = 5
    efficiency: float = 0.8
    bias_sd: float = 0.0
    sampling: str = 'Fixed coverage'
    pcr_sub: float = 0.00001
    pcr_ins: float = 0.0
    pcr_del: float = 0.0
    sequencing: bool = True
    seq_sub: float = 0.003
    seq_ins: float = 0.0002
    seq_del: float = 0.0002
    quality_mode: str = "Error rate"
    q_profile: str = "Fixed Q"
    q_fixed: int = 30
    q_start: int = 35
    q_end: int = 20
    q_std: float = 3.0
    paired: bool = False
    read_length: int = 150
    max_variants: int = 100000
    abundance_sigma: float = 0.0
    storage_survival: float = 1.0
    storage_halflives: float = 0.0
    storage_rate_per_halflife: float = 0.0
    efficiency_relative_sd: float = 0.0


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def mutate(s, rng, sub=0., ins=0., delete=0., ct=0.):
    out, events = [], Counter()
    for b in s:
        if rng.random() < ins:
            out.append(str(rng.choice(list("ACGT")))); events["ins"] += 1
        if rng.random() < delete:
            events["del"] += 1; continue
        if b == "C" and rng.random() < ct:
            b = "T"; events["C_to_T"] += 1; events["sub"] += 1
        if rng.random() < sub:
            b = str(rng.choice([x for x in "ACGT" if x != b])); events["sub"] += 1
        out.append(b)
    return "".join(out), events


def simulate(rows, cfg, progress=None, stage_progress=None, stop_after=None, checkpoint=None):
    # Local import: measurements run only in the simulator, never in decoding.
    from evaluation import edits
    if not 0<=cfg.bias_sd<=.5 or cfg.sampling not in ('Fixed coverage','Pool sampling','Pool sampling + minimum 1'):
        raise ValueError('Invalid PCR bias or sampling mode')
    if cfg.scope not in ("BODY only", "Full strand"):
        raise ValueError("Unknown error scope")
    if not rows or not 1 <= cfg.founders <= 100 or not 1 <= cfg.coverage <= 10000:
        raise ValueError("Require oligos, 1..100 founders and 1..10000 reads per oligo")
    if not 0 <= cfg.cycles <= 20 or not 1 <= cfg.read_length <= 300:
        raise ValueError("PCR cycles or read length outside supported range")
    if len(rows)*cfg.founders > cfg.max_variants:
        raise ValueError("Too many founder molecules. Reduce founders or library size.")
    for name in ("syn_sub", "syn_ins", "syn_del", "storage_sub", "storage_ct", "efficiency", "pcr_sub", "pcr_ins", "pcr_del", "seq_sub", "seq_ins", "seq_del"):
        if not 0 <= getattr(cfg,name) <= 1:
            raise ValueError(f"Invalid probability: {name}")
    if not 0 <= cfg.abundance_sigma <= 3 or not 0 <= cfg.storage_survival <= 1 or not 0 <= cfg.efficiency_relative_sd <= .5:
        raise ValueError("Invalid abundance, survival or relative PCR variation")
    if not 0<=cfg.storage_halflives<=10 or not 0<=cfg.storage_rate_per_halflife<=.1:
        raise ValueError("Invalid storage half-life exposure or damage rate")
    rng = np.random.default_rng(cfg.seed)
    outputs = [[], []]; truth, stages, trace = [], [], []
    stage_totals = {}; snapshots = Counter(); stage_quality=[]; populations=[]; abundance=[]
    audit=[]
    order_names=['Synthesis','Storage','PCR','Sequencing']
    start=0
    if checkpoint:
        checkpoint=copy.deepcopy(checkpoint)
        start=order_names.index(checkpoint['completed_stage'])+1
        rng.bit_generator.state=checkpoint['_rng_state']
        stage_totals=checkpoint['_stage_totals'];trace=checkpoint['trace']
        stage_quality=checkpoint['stage_quality'];abundance=checkpoint['abundance'];audit=checkpoint.get('audit',[])
    def audit_pair(stage,before,after,events,founder,cycle=None):
        # Bounded examples per strand/stage; aggregate event counters remain complete.
        key=(stage,sid)
        if audit_counts[key]>=4:
            if not any(events.get(k,0) for k in ('sub','ins','del')):return
            replace_index=next((i for i,r in enumerate(audit) if (r['Stage'],r['Strand ID'])==key and not (r['New S']+r['New I']+r['New D'])),None)
            if replace_index is None:return
            audit.pop(replace_index)
        else:audit_counts[key]+=1
        audit.append({'Stage':stage,'Strand ID':sid,'Founder':founder,'Cycle':cycle,
                      'Before':before,'After':after,'Before nt':len(before),'After nt':len(after),
                      'New S':events.get('sub',0),'New I':events.get('ins',0),'New D':events.get('del',0)})
    audit_counts=Counter((r['Stage'],r['Strand ID']) for r in audit)
    def save_stage(name):
        return {'completed_stage':name,'_populations':copy.deepcopy(populations),
                '_rng_state':copy.deepcopy(rng.bit_generator.state),'_stage_totals':copy.deepcopy(stage_totals),
                'settings':asdict(cfg),'stages':[{'Stage':n,**dict(c)} for n,c in stage_totals.items()],
                'stage_quality':copy.deepcopy(stage_quality),'trace':copy.deepcopy(trace),
                'audit':copy.deepcopy(audit),'abundance':copy.deepcopy(abundance),
                '_initial_counts':initial_counts.copy()}
    last_progress={}
    def notify(stage,fraction):
        bucket=int(fraction*100)
        if stage_progress and last_progress.get(stage)!=bucket:stage_progress(stage,fraction)
        last_progress[stage]=bucket
    populations=[Counter({(f,str(row['Payload']) if cfg.scope=='BODY only' else str(row['Full strand'])):1 for f in range(cfg.founders)}) for row in rows]
    if checkpoint:
        populations=checkpoint["_populations"]
        initial_counts=checkpoint.get('_initial_counts',[cfg.founders]*len(rows))
    else:
        initial_counts=[cfg.founders]*len(rows)
        if cfg.abundance_sigma:
            # Fixed total, unequal allocation. At least one molecule per SID is
            # explicitly reserved: this is a no-synthesis-dropout approximation.
            weights=rng.lognormal(-cfg.abundance_sigma**2/2,cfg.abundance_sigma,len(rows))
            initial_counts=(1+rng.multinomial(len(rows)*(cfg.founders-1),weights/weights.sum())).tolist()
            populations=[Counter({(f,str(row['Payload']) if cfg.scope=='BODY only' else str(row['Full strand'])):1
                                  for f in range(n)}) for row,n in zip(rows,initial_counts)]
    def snapshot(stage, events):
        total = stage_totals.setdefault(stage, Counter())
        total.update(events); total["molecules"] += sum(variants.values())
        total["variants"] += len(variants); total["oligos"] += 1
        for (founder,seq), n in variants.items():
            if row_index < 3 and founder == 0:
                trace.append({"sid":sid,"founder":founder,"stage":stage,"sequence":seq,"copies":n})
        snapshots[stage] += 1
    for row_index,row in enumerate(rows):
        sid=int(row['Strand ID']);variants=populations[row_index]
        if not checkpoint:snapshot('Designed founders',Counter())
    for stage,enabled,params in (('Synthesis',cfg.synthesis,(cfg.syn_sub,cfg.syn_ins,cfg.syn_del,0.)),('Storage',cfg.storage,(1-(1-cfg.storage_sub)*(1-min(1.,cfg.storage_rate_per_halflife*cfg.storage_halflives)),0.,0.,cfg.storage_ct))):
        if order_names.index(stage)<start:continue
        notify(stage,0.)
        for row_index,row in enumerate(rows):
            sid=int(row['Strand ID']);variants=populations[row_index]
            clean=str(row['Payload']) if cfg.scope=='BODY only' else str(row['Full strand'])
            def measure(population):
                counts=Counter();exact=0
                for (_founder,sequence),copies in population.items():
                    e=edits(clean,sequence)
                    for key in ('S','I','D','errors'):counts[key]+=copies*e[key]
                    exact+=copies*(sequence==clean)
                count=sum(population.values());den=count*len(clean)
                return {**dict(counts),'molecules':count,'reference_nt':den,
                        'error_rate_percent':100*counts['errors']/den if den else None,
                        'exact_molecules_percent':100*exact/count if count else None}
            before=measure(variants)
            new, ev = Counter(), Counter()
            for (founder,seq), n in variants.items():
                if stage=='Storage' and enabled and cfg.storage_survival<1:
                    survived=int(rng.binomial(n,cfg.storage_survival));ev['lost_molecules']+=n-survived;n=survived
                for _ in range(n):
                    changed, e = mutate(seq,rng,*params) if enabled else (seq,Counter())
                    new[(founder,changed)] += 1; ev.update(e)
                    ev['input_nt']+=len(seq)
                    audit_pair(stage,seq,changed,e,founder)
            variants = new; snapshot(stage,ev)
            stage_quality.append({'Stage':stage,'Strand ID':sid,'Enabled':enabled,
                                  **{'Before '+k:v for k,v in before.items()},
                                  **{'After '+k:v for k,v in measure(variants).items()}})
            populations[row_index]=variants
            notify(stage,(row_index+1)/len(rows))
        if stop_after==stage:return save_stage(stage)
    if start<=2:notify('PCR',0.)
    for row_index,row in enumerate(rows if start<=2 else []):
        sid=int(row['Strand ID']);variants=populations[row_index]
        efficiency=float(np.clip(rng.normal(cfg.efficiency,cfg.bias_sd),0,1)) if cfg.pcr and cfg.bias_sd else cfg.efficiency
        if cfg.pcr and cfg.efficiency_relative_sd:
            efficiency=float(np.clip(rng.normal(1.,cfg.efficiency_relative_sd)*(1+cfg.efficiency)-1,0,1))
        for cycle in range(cfg.cycles if cfg.pcr else 0):
            new, ev, work = variants.copy(), Counter(), 0
            for (founder,seq), n in variants.items():
                children = int(rng.binomial(n,efficiency))
                ev["input_nt"]+=children*len(seq)
                # Event independence: sub/ins/del Bernoulli at each current base.
                no_event = ((1-cfg.pcr_sub)*(1-cfg.pcr_ins)*(1-cfg.pcr_del))**len(seq)
                changed = int(rng.binomial(children,1-no_event))
                new[(founder,seq)] += children-changed
                if children>changed:audit_pair("PCR",seq,seq,Counter(),founder,cycle+1)
                work += changed
                if work > cfg.max_variants:
                    raise ValueError("PCR work budget exceeded. Lower cycles/error rate; no silent truncation.")
                for _ in range(changed):
                    # Sample conditional on at least one event, without repeated
                    # rejection loops at very small error probabilities.
                    probabilities = [cfg.pcr_sub,cfg.pcr_ins,cfg.pcr_del]*len(seq)
                    first_weights=[]; survival=1.0
                    for p in probabilities:
                        first_weights.append(survival*p); survival*=1-p
                    weights=np.asarray(first_weights); weights/=weights.sum()
                    first=int(rng.choice(len(weights),p=weights))
                    active=np.zeros(len(weights),dtype=bool); active[first]=True
                    for j in range(first+1,len(weights)):
                        active[j]=rng.random()<probabilities[j]
                    chars=[];child_events=Counter()
                    for j,b in enumerate(seq):
                        sub,ins,dele=active[3*j:3*j+3]
                        if ins:
                            chars.append(str(rng.choice(list("ACGT")))); ev["ins"]+=1;child_events["ins"]+=1
                        if dele:
                            ev["del"]+=1;child_events["del"]+=1; continue
                        if sub:
                            b=str(rng.choice([x for x in "ACGT" if x!=b])); ev["sub"]+=1;child_events["sub"]+=1
                        chars.append(b)
                    new[(founder,"".join(chars))]+=1
                    audit_pair("PCR",seq,"".join(chars),child_events,founder,cycle+1)
            variants=+new
            if len(variants)>cfg.max_variants:
                raise ValueError("PCR variant budget exceeded; reduce settings")
            snapshot(f"PCR cycle {cycle+1}",ev)
        if sum(len(v) for v in populations)-len(populations[row_index])+len(variants)>cfg.max_variants:
            raise ValueError('Total PCR variant budget exceeded; reduce cycles or error rates')
        populations[row_index]=variants
        abundance.append({'Strand ID':sid,'Initial molecules':initial_counts[row_index],'PCR efficiency':efficiency if cfg.pcr else None,'Final molecules':sum(variants.values())})
        notify("PCR",(row_index+1)/len(rows))
        if progress:progress(.5*(row_index+1)/len(rows))
    if stop_after=="PCR":return save_stage("PCR")
    totals=np.array([sum(v.values()) for v in populations],float)
    if not totals.sum():raise ValueError('No surviving molecules to sequence')
    if np.any(totals==0) and cfg.sampling!='Pool sampling':raise ValueError('Cannot guarantee coverage for a missing strand')
    probabilities=totals/totals.sum()
    if cfg.sampling=='Fixed coverage':allocations=np.full(len(rows),cfg.coverage,dtype=int)
    else:
        minimum=1 if cfg.sampling=='Pool sampling + minimum 1' else 0
        allocations=minimum+rng.multinomial(len(rows)*(cfg.coverage-minimum),probabilities)
    notify('Sequencing',0.)
    for row_index,(row,variants,nreads) in enumerate(zip(rows,populations,allocations)):
        sid=int(row['Strand ID']);prefix=str(row['FBR'])+str(row['SI']);suffix=str(row.get('Filler',''))+str(row['RBR'])
        abundance[row_index]['Read units']=int(nreads)
        if not nreads:continue
        keys=list(variants); weights=np.array([variants[k] for k in keys],float); weights/=weights.sum()
        for index in rng.choice(len(keys),int(nreads),p=weights):
            founder,template=keys[int(index)]
            if cfg.scope=="BODY only":
                full=prefix+template+suffix
                susceptible=[False]*len(prefix)+[True]*len(template)+[False]*len(suffix)
            else:
                full=template; susceptible=[True]*len(full)
            # No design SID, founder, true boundary or mutation flag in FASTQ.
            unit=len(truth); truth.append({"unit":unit,"sid":sid,"founder":int(founder),"template":full})
            for mate in range(2 if cfg.paired else 1):
                seq=full if mate==0 else rc(full)
                mask=susceptible if mate==0 else susceptible[::-1]
                seq=seq[:cfg.read_length]; mask=mask[:cfg.read_length]
                out, qs, ev=[],[],Counter()
                ev["input_nt"]=sum(mask)
                ev["read_template_nt"]=len(seq)
                for pos,(b,eligible) in enumerate(zip(seq,mask)):
                    if cfg.quality_mode=="Error rate":
                        q=round(-10*math.log10(max(cfg.seq_sub,1e-6)))
                    elif cfg.q_profile=="Linear quality decay":
                        q=cfg.q_start+(cfg.q_end-cfg.q_start)*pos/max(1,len(seq)-1)
                    elif cfg.q_profile=="Normal Q distribution":
                        q=rng.normal(cfg.q_fixed,cfg.q_std)
                    else:
                        q=cfg.q_fixed
                    q=int(np.clip(round(q),2,60))
                    if not cfg.sequencing: q=60
                    psub=cfg.seq_sub if cfg.quality_mode=="Error rate" else 10**(-q/10)
                    if eligible and cfg.sequencing:
                        if rng.random()<cfg.seq_ins:
                            out.append(str(rng.choice(list("ACGT")))); qs.append(chr(q+33)); ev["ins"]+=1
                        if rng.random()<cfg.seq_del:
                            ev["del"]+=1; continue
                        if rng.random()<psub:
                            b=str(rng.choice([x for x in "ACGT" if x!=b])); ev["sub"]+=1
                    out.append(b); qs.append(chr(q+33))
                audit_pair("Sequencing",seq,"".join(out),ev,founder)
                if not out:
                    raise ValueError("A read became empty under extreme errors. Reduce rates; no read was silently replaced.")
                outputs[mate].append(f'@read{unit:09d}/{mate+1}\n'+"".join(out)+"\n+\n"+"".join(qs)+"\n")
                t=stage_totals.setdefault("Sequencing",Counter()); t.update(ev); t["reads"]+=1
        notify("Sequencing",(row_index+1)/len(rows))
        if progress: progress(.5+.5*(row_index+1)/len(rows))
    for name,count in stage_totals.items():
        stages.append({"Stage":name,**dict(count)})
    order=rng.permutation(len(truth))
    for mate in range(2 if cfg.paired else 1):
        outputs[mate]=[f'@read{i:09d}/{mate+1}\n'+outputs[mate][int(old)].split('\n',1)[1] for i,old in enumerate(order)]
    truth=[{**truth[int(old)],'unit':i} for i,old in enumerate(order)]
    return {"r1":"".join(outputs[0]).encode(),"r2":"".join(outputs[1]).encode(),
            "settings":asdict(cfg),"stages":stages,"stage_quality":stage_quality,"abundance":abundance,"trace":trace,"truth":truth,
            "audit":audit,"completed_stage":"Sequencing",
            "note":"Sampling mode is explicit. Pool sampling can yield zero reads without physical molecule dropout. Defaults are synthetic, not vendor-calibrated."}
