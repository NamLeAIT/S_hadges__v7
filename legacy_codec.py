from __future__ import annotations
from typing import List, Dict, Tuple
import re
import hashlib

BASES = 'ACGT'
DIMERS = [a+b for a in BASES for b in BASES]
_SIMPLE_ENC = {'00':'A','01':'C','10':'G','11':'T'}
_SIMPLE_DEC = {v:k for k,v in _SIMPLE_ENC.items()}

DISPLAY = {
    'Simple Mapping':'SM',
    'R0_B9':'R0-B9 (paper)',
    'R1_B12':'R1-B12 (paper)',
    'R2_B15':'R2-B15 (paper)',
    'RINF_B16':'R∞-B16 (paper)',
    'RINF_SHUFFLED_B16':'R∞',
}
RULE_BASE = {'R0_B9':9,'R1_B12':12,'R2_B15':15,'RINF_B16':16}
RULE_BASE['RINF_SHUFFLED_B16']=16
INIT_DIMER = 'TA'


def clean_dna_text(s: str) -> str:
    return re.sub(r'[^ACGTacgt]', '', str(s or '')).upper()


def gc_content(dna: str) -> float:
    dna = clean_dna_text(dna)
    return sum(c in 'GC' for c in dna)/len(dna) if dna else 0.0


def longest_homopolymer(dna: str) -> int:
    dna = clean_dna_text(dna)
    if not dna: return 0
    best = cur = 1
    for i in range(1, len(dna)):
        if dna[i] == dna[i-1]:
            cur += 1; best = max(best, cur)
        else:
            cur = 1
    return best


def bytes_to_bits(data: bytes) -> str:
    return ''.join(f'{b:08b}' for b in bytes(data or b''))


def bits_to_bytes(bits: str) -> bytes:
    if any(c not in '01' for c in bits):
        raise ValueError('bits must contain only 0/1')
    if len(bits) % 8:
        raise ValueError('decoded bitstream is not byte-aligned')
    return bytes(int(bits[i:i+8],2) for i in range(0,len(bits),8))


def simple_encode_bits_to_dna(bits: str) -> Tuple[str, List[int]]:
    if any(c not in '01' for c in bits):
        raise ValueError('bits must be 0/1')
    pad = len(bits) % 2
    payload = bits + ('0' if pad else '')
    header = '01' if pad else '00'
    chunks = [header] + [payload[i:i+2] for i in range(0,len(payload),2)]
    return ''.join(_SIMPLE_ENC[ch] for ch in chunks), [int(ch,2) for ch in chunks]


def simple_decode_dna_to_bits(dna: str) -> Tuple[str, List[int]]:
    dna = clean_dna_text(dna)
    chunks = [_SIMPLE_DEC[b] for b in dna]
    if not chunks: return '', []
    header = chunks[0]
    if header not in {'00','01'}: raise ValueError('Corrupted SIMPLE stream header')
    bits = ''.join(chunks[1:])
    if header == '01':
        if not bits: raise ValueError('Corrupted SIMPLE padding header')
        bits = bits[:-1]
    return bits, [int(ch,2) for ch in chunks]


def allowed_dimers(name: str, prev: str) -> List[str]:
    if name in ('RINF_B16','RINF_SHUFFLED_B16'):
        return list(DIMERS)
    if name == 'R2_B15':
        y = prev[1]
        return [d for d in DIMERS if d != y+y]
    if name == 'R1_B12':
        return [d for d in DIMERS if d[0] != d[1]]
    if name == 'R0_B9':
        y = prev[1]
        return [a+b for a in BASES if a != y for b in BASES if a != b]
    raise ValueError(f'Unknown rule {name}')


def bits_to_base_digits(bits: str, base: int, prepend_one: bool=True) -> List[int]:
    if not bits or any(c not in '01' for c in bits):
        raise ValueError("bits must be a non-empty string of '0'/'1'")
    s = ('1'+bits) if prepend_one else bits
    n = int(s,2)
    out=[]
    while n>0:
        out.append(n%base); n//=base
    return out[::-1] if out else [0]


def base_digits_to_bits(digits: List[int], base: int, remove_leading_one: bool=True) -> str:
    if not digits: raise ValueError('digits cannot be empty')
    n=0
    for d in digits:
        if not (0 <= d < base): raise ValueError(f'digit {d} out of range for base={base}')
        n=n*base+d
    b=bin(n)[2:]
    if remove_leading_one:
        if not b or b[0] != '1': raise ValueError("Corrupted stream: leading '1' missing")
        b=b[1:]
    return b


def rule_order(name,prev,step):
    if name=='RINF_SHUFFLED_B16':
        return sorted(DIMERS,key=lambda d:hashlib.sha256(f'rinf-shuffle-v1:{step}:{prev}:{d}'.encode('ascii')).digest())
    return sorted(allowed_dimers(name,prev))


def encode_rule_bits(bits: str, name: str, init_dimer: str=INIT_DIMER) -> Tuple[str,List[int]]:
    base=RULE_BASE[name]
    digits=bits_to_base_digits(bits if bits else '0', base, prepend_one=True)
    prev=init_dimer
    out=[init_dimer]
    for step,d in enumerate(digits,1):
        allowed=rule_order(name,prev,step)
        if len(allowed) != base: raise RuntimeError('Rule table cardinality mismatch')
        nxt=allowed[d]
        out.append(nxt); prev=nxt
    return ''.join(out), digits


def decode_rule_dna(dna: str, name: str, init_dimer: str=INIT_DIMER, strict: bool=True) -> Tuple[str,List[int]]:
    dna=clean_dna_text(dna)
    if len(dna)%2: raise ValueError('Rule-based DNA length must be even')
    units=[dna[i:i+2] for i in range(0,len(dna),2)]
    if not units: return '',[]
    observed_init=units[0]
    if strict and observed_init != init_dimer:
        raise ValueError(f'Initiator mismatch: expected {init_dimer}, observed {observed_init}')
    prev=observed_init
    digits=[]
    for pos,current in enumerate(units[1:], start=1):
        allowed=rule_order(name,prev,pos)
        if current not in allowed:
            raise ValueError(f'Invalid dimer {current} at rule step {pos} for previous state {prev}')
        digits.append(allowed.index(current)); prev=current
    bits=base_digits_to_bits(digits,RULE_BASE[name],remove_leading_one=True)
    return bits,digits


def encode_bytes(data: bytes, design: str) -> Tuple[str,str,Dict]:
    bits=bytes_to_bits(data)
    work=bits if bits else '0'
    if design == 'Simple Mapping':
        dna,digits=simple_encode_bits_to_dna(work)
    else:
        dna,digits=encode_rule_bits(work,design)
    return dna,bits,{
        'design':design,'bits_len':len(bits),'bytes_len':len(data),'digits_len':len(digits),
        'family':'simple' if design=='Simple Mapping' else 'paper_rulebase_B_original',
        'error_containment': False if design!='Simple Mapping' else None,
    }


def decode_dna(dna: str, design: str, meta: Dict|None=None) -> Tuple[bytes,str,Dict]:
    meta=meta or {}
    if design == 'Simple Mapping':
        bits,_=simple_decode_dna_to_bits(dna)
    else:
        bits,_=decode_rule_dna(dna,design)
    if 'bits_len' in meta:
        bits=bits[:int(meta['bits_len'])]
    if len(bits)%8:
        raise ValueError('Decoded legacy bitstream is not byte-aligned')
    data=bits_to_bytes(bits)
    if 'bytes_len' in meta:
        data=data[:int(meta['bytes_len'])]
    return data,bits,{'design':design,'decoded_bits_len':len(bits),'decoded_bytes_len':len(data)}
