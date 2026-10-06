"""Literature-derived parameters with explicit approximations, not lab calibration."""
from dataclasses import asdict
import hashlib
from pathlib import Path
import process_model as pm

PAPER='https://doi.org/10.1038/s41467-023-41729-1'
REPO='https://github.com/fml-ethz/dt4dds'
PRESET='dt4dds-derived · material deposition + Q5 (simplified)'


def defaults(seed=12345):
    # Basic mode is deliberately a small-founder teaching/robustness scenario.
    # Advanced keeps the wider 1..100 range for sensitivity analysis.
    return pm.Settings(seed=seed,founders=2,abundance_sigma=.30,scope='Full strand',
        syn_sub=0.,syn_ins=0.,syn_del=.0005695,storage=False,storage_sub=0.,storage_ct=0.,
        storage_survival=1.,storage_rate_per_halflife=.000164,storage_halflives=0.,cycles=10,efficiency=.8,bias_sd=0.,efficiency_relative_sd=.0051,
        pcr_sub=1.09e-4/280.,pcr_ins=0.,pcr_del=0.,sampling='Pool sampling',coverage=350,paired=True,read_length=125)


def provenance():
    return [
      {'Process':'Synthesis','Parameter':'Deletion probability','Value':.0005695,'Unit':'per base / synthesized copy',
       'Source':'array_synthesis_twist.yaml','Status':'Literature-derived approximation',
       'Implementation':'Independent deletion at each base; original base/run biases are NOT reproduced.'},
      {'Process':'Synthesis','Parameter':'Abundance lognormal sigma','Value':.30,'Unit':'natural-log scale',
       'Source':'array_synthesis_twist.yaml','Status':'Literature-derived approximation',
       'Implementation':'Lognormal weights; allocate a fixed molecule budget, reserving one per SID. Not identical to dt4dds.'},
      {'Process':'Synthesis','Parameter':'Substitution / insertion','Value':0.,'Unit':'per base',
       'Source':'array_synthesis_twist.yaml','Status':'Source model choice',
       'Implementation':'Zero is a modeling assumption, not a vendor zero-error guarantee.'},
      {'Process':'PCR','Parameter':'Q5 substitution probability','Value':1.09e-4/280.,'Unit':'per newly copied base',
       'Source':'pcr_Q5.yaml; processes/pcr.py','Status':'Literature-derived approximation',
       'Implementation':'Base rate / relative fidelity. Applied to each child in a branching model; uniform substitution spectrum; not dt4dds aggregate-cycle implementation.'},
      {'Process':'PCR','Parameter':'Relative amplification-factor SD','Value':.0051,'Unit':'relative factor',
       'Source':'properties.yaml; processes/pcr.py','Status':'Literature-derived parameter',
       'Implementation':'e_i=clip((1+e_mean)*Normal(1,SD)-1,0,1), fixed for a SID; inherited by its variants.'},
      {'Process':'Storage','Parameter':'Substitution rate per half-life','Value':.000164,'Unit':'per base per half-life',
       'Source':'aging.yaml; processes/aging.py','Status':'Literature-derived approximation',
       'Implementation':'p=rate * exposure. Uniform substitution spectrum; source chemical/base bias is not reproduced. Optional whole-molecule survival 2**(-exposure) approximates loss, not fragment breakage.'},
      {'Process':'Storage','Parameter':'Default storage exposure','Value':0.,'Unit':'exposure',
       'Source':'User workflow assumption','Status':'Bypassed',
       'Implementation':'Default zero exposure. Exposure is in half-lives, not years. C-to-T override is custom.'},
      {'Process':'Workflow','Parameter':'Founders / PCR cycles / mean efficiency','Value':'2 / 10 / 0.8','Unit':'molecules / cycles / probability',
       'Source':'App scenario choices','Status':'Assumed; user-adjustable',
       'Implementation':'Not claimed as the paper experimental protocol or real synthesis yield.'},
    ]


def manifest(cfg):
    base=asdict(defaults(cfg.seed));current=asdict(cfg)
    relevant=['founders','abundance_sigma','scope','synthesis','syn_sub','syn_ins','syn_del','storage',
              'storage_sub','storage_ct','storage_survival','storage_halflives','storage_rate_per_halflife','pcr','cycles','efficiency','efficiency_relative_sd','pcr_sub','pcr_ins','pcr_del']
    changed={k:{'preset':base[k],'run':current[k]} for k in relevant if base[k]!=current[k]}
    root=Path(__file__).parent/'preset_sources'
    return {'name':PRESET,'label':('Custom — based on ' if changed else '')+PRESET,'paper':PAPER,'repository':REPO,
            'calibration':'Literature-derived simplified composite scenario; NOT a validated dt4dds reproduction or lab calibration.',
            'overrides':changed,'parameters':provenance(),
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob('*.yaml')}}
