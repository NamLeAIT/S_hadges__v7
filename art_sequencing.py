"""Real ART adapter. References/lineage are evaluation-only; FASTQ names are neutral.

No synthetic fallback. Sequence and quality characters from ART are unchanged.
Counts refer to reads (SE) / pairs (PE), never genomic fold coverage.
"""
from dataclasses import dataclass, asdict
from collections import Counter, defaultdict
from pathlib import Path
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
import time
import numpy as np
from process_model import rc
from read_pipeline import parse


@dataclass(frozen=True)
class ArtSettings:
    executable: str = ''
    profile: str = 'HS25'
    read_length: int = 125
    paired: bool = True
    seed: int = 12345
    q_shift_r1: int = 0
    q_shift_r2: int = 0
    short_policy: str = 'Shorten reads to template'
    timeout_seconds: int = 300
    # ART Windows 2.5.8 can emit malformed FASTQ when its own indel simulation
    # changes read length but leaves the original quality-string length.  The
    # default therefore preserves ART's empirical quality/substitution model
    # while disabling ART's sequencing indels. Molecular indels from mode 2
    # are retained in the input template.
    insertion_r1: float | None = 0.0
    deletion_r1: float | None = 0.0
    insertion_r2: float | None = 0.0
    deletion_r2: float | None = 0.0


def detect_art_executable(custom=''):
    value=custom.strip() or os.environ.get('ART_ILLUMINA','').strip()
    if value:
        path=Path(value).expanduser()
        if path.is_file():return str(path.resolve())
        return shutil.which(value)
    return shutil.which('art_illumina') or shutil.which('art_illumina.exe')


def validate(c):
    if c.profile!='HS25':raise ValueError('This release validates the HS25 profile only.')
    if not 25<=c.read_length<=150:raise ValueError('HS25 read length must be 25..150 nt in this release.')
    if c.short_policy not in ('Stop if template is short','Shorten reads to template'):
        raise ValueError('Unknown short-template policy')
    if not 0<=c.seed<2**31 or not 1<=c.timeout_seconds<=3600:
        raise ValueError('Invalid ART seed or timeout')
    if not all(-30<=v<=30 for v in (c.q_shift_r1,c.q_shift_r2)):
        raise ValueError('Quality shifts must be between -30 and +30')
    for v in (c.insertion_r1,c.deletion_r1,c.insertion_r2,c.deletion_r2):
        if v is not None and not 0<=v<=.1:raise ValueError('ART indel probabilities must be 0..0.1')


def pool_templates(rows, checkpoint):
    """Collapse identical sequences per SID, preserving their total abundance."""
    out=[]
    for row,pop in zip(rows,checkpoint['_populations']):
        grouped=Counter()
        for (_founder,seq),count in pop.items():
            full=str(row['FBR'])+str(row['SI'])+seq+str(row.get('Filler',''))+str(row['RBR']) if checkpoint['settings']['scope']=='BODY only' else seq
            if count:grouped[full]+=int(count)
        for seq,count in sorted(grouped.items()):
            out.append({'sid':int(row['Strand ID']),'template':seq,'copies':count})
    return out


def allocate(templates,total,seed,mode='Pool-proportional'):
    """Allocate sequencing units either from the pooled population or equally
    by designed SID.  The latter still samples molecular variants within each
    SID in proportion to their PCR abundance.
    """
    if not templates or sum(x['copies'] for x in templates)<=0:
        raise ValueError('No surviving molecules to sequence.')
    rng=np.random.default_rng(seed)
    if mode=='Pool-proportional':
        weights=np.array([x['copies'] for x in templates],float)
        counts=rng.multinomial(total,weights/weights.sum())
        return [{**x,'count':int(n)} for x,n in zip(templates,counts) if n]
    if mode!='Equal per strand':raise ValueError('Unknown sequencing sampling mode')
    by_sid=defaultdict(list)
    for x in templates:by_sid[x['sid']].append(x)
    sids=sorted(by_sid)
    if not sids:raise ValueError('No surviving strands to sequence.')
    per=np.full(len(sids),total//len(sids),dtype=int)
    per[:total%len(sids)]+=1
    out=[]
    for sid,n in zip(sids,per):
        group=by_sid[sid];weights=np.array([x['copies'] for x in group],float)
        counts=rng.multinomial(int(n),weights/weights.sum())
        out.extend({**x,'count':int(k)} for x,k in zip(group,counts) if k)
    return out


def pair_name(name):
    name=name.split()[0]
    return name[:-2] if name.endswith(('/1','/2')) else name


def sam_records(path):
    result={}
    for line in path.read_text().splitlines():
        if not line or line.startswith('@'):continue
        f=line.split('\t');flag=int(f[1])
        if flag & (4|256|2048):raise ValueError('Unexpected unmapped/secondary ART alignment')
        cigar=re.findall(r'(\d+)([MIDNSHP=X])',f[5])
        if ''.join(n+op for n,op in cigar)!=f[5] or any(op in 'MNSHP' for _,op in cigar):
            raise ValueError('ART must emit explicit =/X/I/D CIGAR for error accounting')
        count=Counter()
        for n,op in cigar:count[op]+=int(n)
        mate=1 if flag & 128 else 0
        key=(pair_name(f[0]),mate)
        if key in result:raise ValueError('Duplicate ART SAM read identifier')
        result[key]={'variant':f[2],'start':int(f[3])-1,'span':count['=']+count['X']+count['D'],
                     'reverse':bool(flag&16),'S':count['X'],'I':count['I'],'D':count['D']}
    return result


def discover_outputs(directory,paired):
    files=sorted([p for p in directory.iterdir() if p.suffix.lower() in ('.fq','.fastq')])
    if len(files)!=(2 if paired else 1):raise ValueError('ART did not produce the expected FASTQ file count')
    if paired:
        first=[p for p in files if re.search(r'(?:1|R1)\.(?:fq|fastq)$',p.name,re.I)]
        second=[p for p in files if re.search(r'(?:2|R2)\.(?:fq|fastq)$',p.name,re.I)]
        if len(first)!=1 or len(second)!=1:raise ValueError('Cannot identify ART R1/R2 outputs')
        files=first+second
    sams=list(directory.glob('*.sam'))
    if len(sams)!=1:raise ValueError('Missing or ambiguous ART SAM evaluation output')
    return files,sams[0]


def run_art(templates,c=ArtSettings(),progress=None):
    """templates contain sid, template, count. No truth enters a decoder."""
    validate(c);exe=detect_art_executable(c.executable)
    if not exe:raise FileNotFoundError('ART executable not found. Install ART and select art_illumina / art_illumina.exe, or set ART_ILLUMINA. No FASTQ was fabricated.')
    total=sum(x['count'] for x in templates)
    if not 1<=total<=100000:raise ValueError('Request 1..100,000 read units per run; split larger experiments.')
    # ART accepts one read length per invocation.  Molecular indels can make
    # sampled templates heterogeneous, so run one ART batch per input length
    # (and per requested count).  Every input template is retained.
    groups=defaultdict(list);lookup={};shortened=0
    for i,x in enumerate(templates):
        if x['count']<=0:continue
        seq=x['template'].upper()
        if not seq or set(seq)-set('ACGT'):raise ValueError('ART input templates must be nonempty A/C/G/T sequences')
        length=c.read_length
        if len(seq)<length:
            if c.short_policy=='Stop if template is short':
                raise ValueError(f'Template length {len(seq)} < requested {length} nt. Choose the explicit shortened-read approximation or lower read length. No templates were removed.')
            length=len(seq);shortened+=x['count']
        if length<25:raise ValueError('A sampled template is shorter than 25 nt; this ART configuration does not support it.')
        ident=f'v{i:08d}';lookup[ident]={**x,'template':seq}
        groups[(int(x['count']),length)].append(ident)
    if len(groups)>2000:raise ValueError('Too many ART batches; reduce coverage/variant diversity.')
    units=[];commands=[];logs=[];stats=Counter();seen=Counter();deadline=time.monotonic()+c.timeout_seconds
    seed_rng=np.random.default_rng(c.seed)
    try:
        with tempfile.TemporaryDirectory(prefix='dna_art_') as tmp:
            for job,((count,length),ids) in enumerate(sorted(groups.items())):
                directory=Path(tmp)/f'batch_{job:04d}';directory.mkdir()
                fasta=directory/'input.fa'
                fasta.write_text(''.join(f'>{i}\n{lookup[i]["template"]}\n' for i in ids),encoding='ascii')
                cmd=[exe,'-ss',c.profile,'-amp','-sam','-na','-i','input.fa','-l',str(length),'-c',str(count),
                     '-rs',str(int(seed_rng.integers(1,2**31-1))),'-qs',str(c.q_shift_r1),'-o','reads_']
                if c.paired:cmd+=['-p','-qs2',str(c.q_shift_r2)]
                for flag,v in [('-ir',c.insertion_r1),('-dr',c.deletion_r1),('-ir2',c.insertion_r2),('-dr2',c.deletion_r2)]:
                    if v is not None:cmd += [flag,str(v)]
                remaining=deadline-time.monotonic()
                if remaining<=0:raise TimeoutError('ART total time budget exceeded; no partial run accepted')
                proc=subprocess.run(cmd,cwd=directory,capture_output=True,text=True,timeout=remaining,check=False)
                logs.append({'batch':job,'stdout':proc.stdout,'stderr':proc.stderr})
                commands.append({'argv':cmd,'template_ids':ids,'input_sha256':hashlib.sha256(fasta.read_bytes()).hexdigest()})
                if proc.returncode:raise RuntimeError(f'ART failed (exit {proc.returncode}): {(proc.stderr+proc.stdout)[-2500:]}')
                files,sam=discover_outputs(directory,c.paired);mates=[parse(p.read_bytes()) for p in files]
                if any(len(m)!=count*len(ids) for m in mates):raise ValueError('ART read count mismatch: a template may have been skipped. No partial result accepted.')
                aligned=sam_records(sam)
                for j,r in enumerate(mates[0]):
                    name=pair_name(r['name']);records=[m[j] for m in mates]
                    if any(pair_name(m['name'])!=name for m in records):raise ValueError('ART mate names/order do not match')
                    spans=[];errors=[];ident=None
                    for mate,rec in enumerate(records):
                        info=aligned.get((name,mate))
                        if not info:raise ValueError('ART FASTQ has no matching SAM alignment')
                        if ident is not None and ident!=info['variant']:raise ValueError('ART pair maps to different templates')
                        ident=info['variant'];template=lookup[ident]['template']
                        span=template[info['start']:info['start']+info['span']]
                        if len(span)!=info['span']:raise ValueError('ART alignment extends past template')
                        spans.append(rc(span) if info['reverse'] else span)
                        errors.append({k:info[k] for k in ('S','I','D')})
                        stats.update({'sub':info['S'],'ins':info['I'],'del':info['D'],'input_nt':info['span'],'reads':1})
                    seen[ident]+=1
                    units.append((records,{'sid':lookup[ident]['sid'],'template':lookup[ident]['template'],
                                          'variant':ident,'sequencing_reference':spans,'sequencing_edits':errors}))
                if progress:progress((job+1)/len(groups))
        if any(seen[i]!=x['count'] for i,x in lookup.items()):raise ValueError('ART per-template count mismatch')
    except subprocess.TimeoutExpired as exc:raise TimeoutError('ART timed out. Increase the Advanced time budget; no partial FASTQ accepted.') from exc
    except OSError as exc:raise RuntimeError(f'Cannot run ART: {exc}. Check executable and runtime dependencies.') from exc
    output=[[],[]];truth=[]
    for i,old in enumerate(seed_rng.permutation(len(units))):
        records,t=units[int(old)];truth.append({**t,'unit':i})
        for mate,r in enumerate(records):
            output[mate].append(f'@read{i:09d}/{mate+1}\n{r["seq"]}\n+\n'+''.join(chr(q+33) for q in r['q'])+'\n')
    length_groups=Counter(len(x['template']) for x in lookup.values())
    metadata={'engine':'ART','settings':asdict(c),'commands':commands,'logs':logs,
              'executable_sha256':hashlib.sha256(Path(exe).read_bytes()).hexdigest(),
              'templates':lookup,'shortened_read_units':shortened,'total_read_units':total,
              'input_template_lengths':dict(sorted(length_groups.items())),
              'short_template_note':'Variable-length input templates were grouped and run in separate ART batches. Shortened reads are an explicit approximation, NOT an adapter/read-through simulation.',
              'pairing_note':'One PE pair represents one sampled template; mates are not independent founder molecules.'}
    return {'r1':''.join(output[0]).encode(),'r2':''.join(output[1]).encode(),'truth':truth,
            'sequencing_stage':{'Stage':'Sequencing',**dict(stats)},'art':metadata}


def run_art_illumina(rows,art_exe='',read_length=125,reads_per_oligo=350,seed=12345,paired=True,**options):
    """Public entry point for mode 1; all upstream molecules are assumed clean."""
    if not 1<=reads_per_oligo<=10000:raise ValueError('Reads/pairs per oligo must be 1..10000')
    templates=[{'sid':int(r['Strand ID']),'template':r['Full strand'],'count':reads_per_oligo} for r in rows]
    return run_art(templates,ArtSettings(executable=art_exe,read_length=read_length,paired=paired,seed=seed,**options))
