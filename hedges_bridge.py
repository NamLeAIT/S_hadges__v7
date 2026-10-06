from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import lru_cache
from typing import Dict, List, Tuple
import random

import hedges_v52 as hv
from framing import bytes_to_bits, frame_103, autonomous_reconstruct, pack_hedges_file, hedges_file_overhead

FBR = "ACACGACGCTCTTCCGATCT"       # accepted App1 primer, 20 nt
RBR = "AGATCGGAAGAGCACACGTCT"      # accepted App1 primer, 21 nt
SI_NT = 8
BODY_NT = 76
FULL_NT = len(FBR) + SI_NT + BODY_NT + len(RBR)
PROFILE_NAME = "P103_T8_C4"
PROFILE = hv.PROFILES[PROFILE_NAME]
MAX_STRANDS = 512
MAX_UNPROTECTED_STRANDS = 8192


def strand_capacity(sid_mode):
    return MAX_UNPROTECTED_STRANDS if normalize_sid_mode(sid_mode)==SID_MODE_NO_HAMMING else MAX_STRANDS
SI_SEED = 20260902
SI_BIO_ONLY_SEED = 20260910
SID_MODE_HAMMING = "Hamming (HD≥3)"
SID_MODE_NO_HAMMING = "No Hamming"
SID_MODE_AUTO = "Auto"
SID_MODES = (SID_MODE_HAMMING, SID_MODE_NO_HAMMING)

assert FULL_NT == 125


@dataclass
class HedgesStrand:
    strand_id: int
    si: str
    payload_bits: str
    body: str
    full_strand: str
    seed_id: int
    encode_status: str
    certified: bool
    rescue: bool
    body_gc: float
    body_hp: int
    full_gc: float
    full_hp: int

    def row(self) -> dict:
        return {
            "Strand ID": self.strand_id,
            "SI": self.si,
            "Payload": self.body,
            "Payload length": BODY_NT,
            "Filler": "",
            "FBR": FBR,
            "RBR": RBR,
            "Full strand": self.full_strand,
            "Full strand length": len(self.full_strand),
            "GC content": self.full_gc,
            "Longest homopolymer": self.full_hp,
            "HEDges seed": self.seed_id,
            "HEDges status": self.encode_status,
            "Certified": self.certified,
            "Rescue": self.rescue,
        }


@dataclass
class HedgesDecoded:
    strand_id: int
    status: str
    payload_bits: str | None
    edit_distance: int | None
    candidate_count: int
    states: int
    selected_seed: int | None


def _hp(seq: str) -> int:
    return int(hv.bio_stats(seq)["longest_homopolymer"])


def _gc(seq: str) -> float:
    return float(hv.bio_stats(seq)["gc_percent"]) / 100.0


def normalize_sid_mode(mode: str | None) -> str:
    value = str(mode or SID_MODE_HAMMING).strip()
    aliases = {
        "Hamming": SID_MODE_HAMMING,
        "Hamming-protected": SID_MODE_HAMMING,
        "HD3": SID_MODE_HAMMING,
        "Bio-constrained only": SID_MODE_NO_HAMMING,
        "No Hamming (bio-constrained)": SID_MODE_NO_HAMMING,
        "No Hamming": SID_MODE_NO_HAMMING,
        "Auto": SID_MODE_AUTO,
    }
    return aliases.get(value, value)


def _si_bio_ok(seq: str) -> bool:
    bs = hv.bio_stats(seq)
    # Match the frozen V5.2 SI-generation biology envelope.  For an 8-nt SI,
    # GC values are discrete, so the historical 25-75% gate is retained.
    return int(bs["longest_homopolymer"]) <= 3 and 25.0 <= float(bs["gc_percent"]) <= 75.0


@lru_cache(maxsize=8)
def _hamming_codebook(n: int) -> Dict[int, str]:
    return hv.generate_si_codebook(n_ids=int(n), min_hd=3, seed=SI_SEED)


@lru_cache(maxsize=8)
def _bio_only_codebook(n: int) -> Dict[int, str]:
    """Deterministic SI codebook with biology constraints but no pairwise HD rule.

    The codebook is kept sequence-disjoint from the Hamming codebook only so an
    uploaded clean SI can identify the public codec mode without hidden session
    state.  No minimum pairwise Hamming distance is imposed among bio-only SIs.
    """
    n = int(n)
    protected = set(_hamming_codebook(MAX_STRANDS).values())
    rng = random.Random(SI_BIO_ONLY_SEED)
    chosen: list[str] = []
    seen: set[str] = set()
    attempts = 0
    while len(chosen) < n and attempts < 2_000_000:
        attempts += 1
        seq = "".join(rng.choice("ACGT") for _ in range(SI_NT))
        if seq in seen or seq in protected:
            continue
        seen.add(seq)
        if not _si_bio_ok(seq):
            continue
        chosen.append(seq)
    if len(chosen) < n:
        raise RuntimeError(f"Could create only {len(chosen)} bio-only SI codewords")
    return {i: seq for i, seq in enumerate(chosen)}


def codebook(n: int = MAX_STRANDS, sid_mode: str = SID_MODE_HAMMING) -> Dict[int, str]:
    mode = normalize_sid_mode(sid_mode)
    if not 1<=n<=strand_capacity(mode):raise ValueError("SI codebook capacity exceeded")
    if mode == SID_MODE_HAMMING:
        return _hamming_codebook(int(n))
    if mode == SID_MODE_NO_HAMMING:
        return _bio_only_codebook(int(n))
    raise ValueError(f"Unsupported SID mode: {sid_mode}")


def decode_si(observed: str, sid_mode: str = SID_MODE_HAMMING) -> Tuple[int, str, int]:
    observed = "".join(c for c in str(observed).upper() if c in "ACGT")
    if len(observed) != SI_NT:
        raise ValueError("SI must be exactly 8 nt")
    mode = normalize_sid_mode(sid_mode)
    if mode == SID_MODE_AUTO:
        return decode_si_auto(observed)[:3]
    cb = codebook(strand_capacity(mode), mode)
    if mode == SID_MODE_NO_HAMMING:
        reverse = {seq: sid for sid, seq in cb.items()}
        if observed not in reverse:
            raise ValueError("SI is not an exact member of the No-Hamming bio-constrained codebook")
        sid = int(reverse[observed])
        return sid, observed, 0

    scored = sorted((hv.hamming(observed, seq), sid, seq) for sid, seq in cb.items())
    d = scored[0][0]
    best = [r for r in scored if r[0] == d]
    if d > 1 or len(best) != 1:
        raise ValueError("SI is not uniquely correctable within HD <= 1")
    _, sid, seq = best[0]
    return int(sid), seq, int(d)


def decode_si_auto(observed: str) -> Tuple[int, str, int, str]:
    """Fail-closed public-mode autodetection for an 8-nt observed SI.

    Exact codeword membership is checked first in both disjoint codebooks. If
    neither matches exactly, only the Hamming-protected codebook is allowed to
    perform its documented one-substitution correction. No-Hamming mode never
    receives nearest-neighbour correction.
    """
    observed = "".join(c for c in str(observed).upper() if c in "ACGT")
    if len(observed) != SI_NT:
        raise ValueError("SI must be exactly 8 nt")
    exact = []
    for mode in SID_MODES:
        cb = codebook(strand_capacity(mode), mode)
        rev = {seq: sid for sid, seq in cb.items()}
        if observed in rev:
            exact.append((int(rev[observed]), observed, 0, mode))
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        raise ValueError("SI mode is ambiguous across public codebooks")
    sid, seq, d = decode_si(observed, SID_MODE_HAMMING)
    return sid, seq, d, SID_MODE_HAMMING


def sid_codebook_stats(n: int, sid_mode: str) -> dict:
    cb = codebook(int(n), sid_mode)
    vals = list(cb.values())
    min_hd = hv.min_pairwise_hamming(cb) if len(vals) > 1 and sid_mode==SID_MODE_HAMMING else None
    gcs = [float(hv.bio_stats(x)["gc_percent"]) for x in vals]
    hps = [int(hv.bio_stats(x)["longest_homopolymer"]) for x in vals]
    return {
        "mode": normalize_sid_mode(sid_mode),
        "count": len(vals),
        "min_pairwise_hamming": min_hd,
        "min_gc_percent": min(gcs) if gcs else None,
        "max_gc_percent": max(gcs) if gcs else None,
        "max_homopolymer": max(hps) if hps else None,
    }


def _biology_key(enc: hv.EncodeResult):
    body = enc.body_bio
    central = enc.central_bio
    outside = max(0.0, 40.0 - body["gc_percent"], body["gc_percent"] - 60.0)
    outside += max(0.0, 40.0 - central["gc_percent"], central["gc_percent"] - 60.0)
    hp = max(int(body["longest_homopolymer"]), int(central["longest_homopolymer"]))
    center = abs(float(body["gc_percent"]) - 50.0) + abs(float(central["gc_percent"]) - 50.0)
    return (outside, hp, center, -hv._selection_tie(enc.payload_bits, enc.salt, PROFILE_NAME, enc.seed_id))


@lru_cache(maxsize=100000)
def select_encoding(payload: str, sid: int, si: str):
    """Deterministic V5.2 selection used by both encoder and decoder.

    Priority: strict certified -> strict mapping-complete feasible -> deterministic rescue.
    """
    strict = hv.select_verified_encoding(
        payload, int(sid), si, PROFILE_NAME,
        seed_count=hv.DEFAULT_SEED_COUNT,
        cert_radius=hv.DEFAULT_CERT_RADIUS,
        strict=True,
    )
    if strict.encoding is not None:
        return strict.encoding, "TVH-V5.2 certified", True, False

    fast = hv.select_fast_encoding(payload, int(sid), si, PROFILE_NAME, hv.DEFAULT_SEED_COUNT)
    if fast.encoding is not None:
        return fast.encoding, "TVH-V5.2 feasible", False, False

    rescue = [
        hv.encode_raw_rescue(payload, int(sid), si, PROFILE, seed)
        for seed in range(hv.RESCUE_SEED_START, hv.RAW_DECODE_HYPOTHESES)
    ]
    enc = min(rescue, key=_biology_key)
    return enc, "TVH-V5.2 rescue", False, True


def estimate_strands(data: bytes, filename: str | None = None, mime: str | None = None) -> int:
    framed_data = pack_hedges_file(data, filename, mime) if filename is not None else bytes(data or b"")
    return (len(framed_data) * 8 + 102) // 103


def encode_bytes(
    data: bytes,
    sid_mode: str = SID_MODE_HAMMING,
    filename: str | None = None,
    mime: str | None = None,
) -> tuple[List[HedgesStrand], dict]:
    source = bytes(data or b"")
    framed_data = pack_hedges_file(source, filename, mime) if filename is not None else source
    bits = bytes_to_bits(framed_data)
    chunks, padding_bits, framed = frame_103(bits)
    limit=strand_capacity(sid_mode)
    if len(chunks) > limit:
        raise ValueError(
            f"TVH-V5.2 supports at most {limit} strands in this SI mode "
            f"({limit * 103 // 8:,} whole bytes approximately)."
        )
    sid_mode = normalize_sid_mode(sid_mode)
    if sid_mode not in SID_MODES:
        raise ValueError(f"Unsupported SID mode: {sid_mode}")
    cb = codebook(strand_capacity(sid_mode), sid_mode)
    out: List[HedgesStrand] = []
    for sid, payload in enumerate(chunks):
        si = cb[sid]
        enc, status, cert, rescue = select_encoding(payload, sid, si)
        full = FBR + si + enc.dna + RBR
        out.append(HedgesStrand(
            strand_id=sid, si=si, payload_bits=payload, body=enc.dna,
            full_strand=full, seed_id=int(enc.seed_id), encode_status=status,
            certified=bool(cert), rescue=bool(rescue),
            body_gc=_gc(enc.dna), body_hp=_hp(enc.dna),
            full_gc=_gc(full), full_hp=_hp(full),
        ))
    return out, {
        "design": "TVH-V5.2 HEDges",
        "sid_mode": sid_mode,
        "sid_codebook": sid_codebook_stats(len(out), sid_mode),
        "payload_bits_per_strand": 103,
        "body_nt": 76,
        "full_strand_nt": 125,
        "strand_count": len(out),
        "padding_bits": padding_bits,
        "framing": "embedded file length/name/domain/CRC32 + self-delimiting tail" if filename is not None else "legacy self-delimiting tail: data|1|0* when padding > 0",
        "original_bits": len(bits),
        "framed_bits": len(framed),
        "original_bytes": len(source),
        "envelope_bytes": len(framed_data),
        "envelope_overhead_bytes": hedges_file_overhead(filename) if filename is not None else 0,
        "stored_filename": filename,
        "stored_mime": mime,
        "rescue_strands": sum(s.rescue for s in out),
        "certified_strands": sum(s.certified for s in out),
    }


def _unique_selected_candidate(raw: hv.DecodeResult, sid: int, si: str) -> HedgesDecoded:
    accepted = {}
    for c in raw.candidates:
        try:
            selected, _status, _cert, _rescue = select_encoding(c.payload_bits, sid, si)
        except Exception:
            continue
        if int(selected.seed_id) != int(c.seed_id):
            continue
        old = accepted.get(c.payload_bits)
        if old is None or (c.edit_distance, c.seed_id) < (old.edit_distance, old.seed_id):
            accepted[c.payload_bits] = c
    vals = sorted(accepted.values(), key=lambda c: (c.edit_distance, c.seed_id, c.payload_bits))
    if len(vals) == 1:
        c = vals[0]
        return HedgesDecoded(sid, "OK", c.payload_bits, int(c.edit_distance), 1, int(raw.states_expanded), int(c.seed_id))
    if not vals:
        return HedgesDecoded(sid, "FAIL", None, None, 0, int(raw.states_expanded), None)
    return HedgesDecoded(sid, "AMBIGUOUS", None, None, len(vals), int(raw.states_expanded), None)


def decode_body(body: str, sid: int, max_edits: int = 2, sid_mode: str = SID_MODE_HAMMING) -> HedgesDecoded:
    sid_mode = normalize_sid_mode(sid_mode)
    cb = codebook(strand_capacity(sid_mode), sid_mode)
    if sid not in cb:
        return HedgesDecoded(sid, "FAIL", None, None, 0, 0, None)
    si = cb[sid]
    raw = hv.raw_all_modes_decode(body, sid, si, PROFILE, int(max_edits))
    return _unique_selected_candidate(raw, sid, si)


def decode_body_auto(body: str, sid: int, max_edits: int = 2) -> HedgesDecoded:
    """Try both public SID codebooks and accept only one unique payload."""
    successful = []
    total_states = 0
    for mode in SID_MODES:
        r = decode_body(body, sid, max_edits=max_edits, sid_mode=mode)
        total_states += int(r.states)
        if r.status == "OK" and r.payload_bits is not None:
            successful.append((mode, r))
    by_payload = {}
    for mode, r in successful:
        by_payload.setdefault(r.payload_bits, []).append((mode, r))
    if len(by_payload) == 1:
        items = next(iter(by_payload.values()))
        best = min(items, key=lambda x: (x[1].edit_distance if x[1].edit_distance is not None else 10**9, x[0]))[1]
        return HedgesDecoded(sid, "OK", best.payload_bits, best.edit_distance, 1, total_states, best.selected_seed)
    if not by_payload:
        return HedgesDecoded(sid, "FAIL", None, None, 0, total_states, None)
    return HedgesDecoded(sid, "AMBIGUOUS", None, None, len(by_payload), total_states, None)


def decode_bodies(bodies: Dict[int, str], max_edits: int = 2, sid_mode: str = SID_MODE_HAMMING):
    mode = normalize_sid_mode(sid_mode)
    results: List[HedgesDecoded] = []
    for sid in sorted(bodies):
        if mode == SID_MODE_AUTO:
            results.append(decode_body_auto(bodies[sid], sid, max_edits=max_edits))
        else:
            results.append(decode_body(bodies[sid], sid, max_edits=max_edits, sid_mode=mode))
    full_bits = "".join(r.payload_bits if r.payload_bits is not None else "?" * 103 for r in results)
    if "?" in full_bits:
        return results, full_bits, "FAIL", None, []
    status, chosen, candidates = autonomous_reconstruct(full_bits)
    return results, full_bits, status, chosen, candidates


def parse_raw_full_strands(text: str, sid_mode: str = SID_MODE_AUTO) -> tuple[Dict[int, str], List[dict]]:
    lines = []
    for raw in str(text).splitlines():
        line = "".join(c for c in raw.upper() if c in "ACGT")
        if line:
            lines.append(line)
    if not lines:
        raise ValueError("No DNA strands found")
    bodies: Dict[int, str] = {}
    rows = []
    for idx, full in enumerate(lines, 1):
        if not full.startswith(FBR) or not full.endswith(RBR):
            raise ValueError(f"Line {idx}: FBR/RBR boundary mismatch")
        central = full[len(FBR):-len(RBR)]
        if len(central) <= SI_NT:
            raise ValueError(f"Line {idx}: strand is too short to contain SI + BODY")
        obs_si = central[:SI_NT]
        mode = normalize_sid_mode(sid_mode)
        if mode == SID_MODE_AUTO:
            sid, canonical, d, detected_mode = decode_si_auto(obs_si)
        else:
            sid, canonical, d = decode_si(obs_si, mode)
            detected_mode = mode
        if sid in bodies:
            raise ValueError(f"Duplicate strand ID {sid}")
        body = central[SI_NT:]
        bodies[sid] = body
        rows.append({"Input line": idx, "Strand ID": sid, "Observed SI": obs_si, "SI distance": d, "SID mode": detected_mode, "BODY nt": len(body)})
    if sorted(bodies) != list(range(len(bodies))):
        raise ValueError("Recovered strand IDs are not contiguous from 0")
    return bodies, rows


def body_reads_for_fastq(strands: List[HedgesStrand]) -> List[dict]:
    """Legacy compatibility helper; V4 UI does not use BODY-only FASTQ simulation.

    V4 simulates full 125-nt strands and resolves read groups from embedded SI.
    """
    return [{"Strand ID": s.strand_id, "Full strand": s.body} for s in strands]
