"""
Terminal-Verified HEDGES V5.2 PARALLEL (TVH-V5.2)

Frozen candidate after V4 stress testing.

BODY = 76 nt = 68-nt main HEDGES region + 8-nt terminal verifier tail.
Useful payload = 103 bits. Internal checkpoints = 4 bits (2+2).
Default seed search = 16; clean-codeword certification radius = 1.
Useful BODY density = 103/76 = 1.355263 bit/nt.
Nominal full-125 density = 103/125 = 0.824.

The 8 terminal bases carry no new payload bits. They are deterministic,
payload/ID/seed/state-dependent verification symbols that supply future
evidence after the last source-bearing position. Decoder is fail-closed:
0 candidates => FAIL, >1 => AMBIGUOUS; never nearest-best unique guess.
V5.2 preserves V5/V5.1 encoding/codewords exactly and changes only decoder execution.

Decoder optimization:
Once all source bits are known, every remaining main/terminal nucleotide is
deterministic. V5.1 generates that entire suffix and closes the path with one
bounded Levenshtein calculation instead of branching through S/I/D hypotheses
one terminal symbol at a time. This is exact, not beam/heuristic pruning.

This is an experimental codec, not a formal arbitrary-K IDS guarantee.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import blake2s
from collections import defaultdict
from functools import lru_cache
import math
import random
import time
import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict

BASES = "ACGT"
SI_NT = 8
BODY_NT = 76
TERMINAL_NT = 8
MAIN_NT = BODY_NT - TERMINAL_NT
HISTORY_BITS = 20
DNA_SUFFIX_NT = 4
ANCHOR_POSITIONS = (6, 15, 25, 36, 48, 59, 65)
DEFAULT_SEED_COUNT = 16
DEFAULT_CERT_RADIUS = 1

MAP_DOMAIN = b"TVH-V5|MAP"  # unchanged: codewords identical to V5
CHECK_DOMAIN = b"TVH-V5|CHECK"
TIE_DOMAIN = b"TVH-V5|TIE"
SELECT_DOMAIN = b"TVH-V5|SELECT"
TERMINAL_DOMAIN = b"TVH-V5|TERMINAL"


@dataclass(frozen=True)
class Profile:
    name: str
    payload_bits: int
    check_widths: tuple[int, ...]

    @property
    def check_bits(self) -> int:
        return sum(self.check_widths)

    @property
    def source_bits(self) -> int:
        return self.payload_bits + self.check_bits

    @property
    def body_density(self) -> float:
        return self.payload_bits / BODY_NT

    @property
    def full125_density(self) -> float:
        return self.payload_bits / 125.0


PROFILES = {
    "P103_T8_C4": Profile("P103_T8_C4", 103, (2, 2)),
}


@dataclass(frozen=True)
class BioState:
    dna_suffix: str = ""
    last_base: str = ""
    run_len: int = 0
    gc_count: int = 0


def _update_bio(state: BioState, base: str) -> BioState:
    run = state.run_len + 1 if base == state.last_base else 1
    return BioState(
        dna_suffix=(state.dna_suffix + base)[-DNA_SUFFIX_NT:],
        last_base=base,
        run_len=run,
        gc_count=state.gc_count + int(base in "GC"),
    )


def initial_bio_from_si(si_seq: str) -> BioState:
    last = ""
    run = 0
    suffix = ""
    for b in si_seq[-DNA_SUFFIX_NT:]:
        run = run + 1 if b == last else 1
        last = b
        suffix = (suffix + b)[-DNA_SUFFIX_NT:]
    return BioState(suffix, last, run, 0)


def bio_stats(dna: str) -> dict:
    if not dna:
        return {"gc_percent": 0.0, "longest_homopolymer": 0, "overall_bio_pass": False}
    gc = 100.0 * sum(b in "GC" for b in dna) / len(dna)
    longest = 0
    last = ""
    run = 0
    for b in dna:
        run = run + 1 if b == last else 1
        last = b
        longest = max(longest, run)
    return {
        "gc_percent": gc,
        "longest_homopolymer": longest,
        "overall_bio_pass": (40.0 <= gc <= 60.0 and longest <= 3),
    }


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal length")
    return sum(x != y for x, y in zip(a, b))


def levenshtein(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j-1] + (ca != cb)))
        prev = cur
    return prev[-1]


def generate_si_codebook(n_ids=256, min_hd=3, seed=20260902, max_attempts=2_000_000):
    if n_ids > 4**SI_NT:
        raise ValueError("8 nt cannot encode >4^8 distinct sequences")
    rng = random.Random(seed)
    chosen = []
    seen = set()
    attempts = 0
    while len(chosen) < n_ids and attempts < max_attempts:
        attempts += 1
        s = "".join(rng.choice(BASES) for _ in range(SI_NT))
        if s in seen:
            continue
        seen.add(s)
        bs = bio_stats(s)
        if bs["longest_homopolymer"] > 3 or not (25.0 <= bs["gc_percent"] <= 75.0):
            continue
        if all(hamming(s, old) >= min_hd for old in chosen):
            chosen.append(s)
    if len(chosen) < n_ids:
        raise RuntimeError(f"Could create only {len(chosen)} SI codewords")
    return {i: s for i, s in enumerate(chosen)}


def min_pairwise_hamming(codebook: Dict[int, str]) -> int | None:
    vals = list(codebook.values())
    if len(vals) < 2:
        return None
    return min(hamming(vals[i], vals[j]) for i in range(len(vals)) for j in range(i+1, len(vals)))


def bits_to_bytes(bits: str) -> bytes:
    if not bits:
        return b""
    return int(bits, 2).to_bytes((len(bits)+7)//8, "big")


def _chunk_lengths(payload_bits: int, n_chunks: int) -> tuple[int, ...]:
    q, r = divmod(payload_bits, n_chunks)
    return tuple(q + int(i < r) for i in range(n_chunks))


@lru_cache(maxsize=500000)
def checkpoint_bits(cumulative_payload: str, stage: int, width: int, salt: int, profile_name: str) -> str:
    material = (
        CHECK_DOMAIN + profile_name.encode() + int(salt).to_bytes(8, "big")
        + int(stage).to_bytes(1, "big") + len(cumulative_payload).to_bytes(2, "big")
        + bits_to_bytes(cumulative_payload)
    )
    digest = blake2s(material, digest_size=4).digest()
    value = int.from_bytes(digest, "big") >> (32-width)
    return f"{value:0{width}b}"


def build_source(payload: str, salt: int, profile: Profile) -> str:
    if len(payload) != profile.payload_bits or any(c not in "01" for c in payload):
        raise ValueError(f"Expected {profile.payload_bits} payload bits")
    lens = _chunk_lengths(profile.payload_bits, len(profile.check_widths))
    parts = []
    cursor = 0
    cumulative = ""
    for stage, (n, width) in enumerate(zip(lens, profile.check_widths)):
        seg = payload[cursor:cursor+n]
        cursor += n
        cumulative += seg
        parts.append(seg)
        parts.append(checkpoint_bits(cumulative, stage, width, salt, profile.name))
    source = "".join(parts)
    assert len(source) == profile.source_bits
    return source


def parse_and_validate_prefix(source_prefix: str, salt: int, profile: Profile):
    if any(c not in "01" for c in source_prefix):
        return False, ""
    lens = _chunk_lengths(profile.payload_bits, len(profile.check_widths))
    cursor = 0
    payload = ""
    for stage, (n, width) in enumerate(zip(lens, profile.check_widths)):
        avail = min(n, max(0, len(source_prefix)-cursor))
        if avail:
            payload += source_prefix[cursor:cursor+avail]
        cursor += avail
        if avail < n:
            return True, payload
        check_avail = min(width, max(0, len(source_prefix)-cursor))
        if check_avail:
            expected = checkpoint_bits(payload, stage, width, salt, profile.name)
            observed = source_prefix[cursor:cursor+check_avail]
            if observed != expected[:check_avail]:
                return False, payload
        cursor += check_avail
        if check_avail < width:
            return True, payload
    return True, payload


@lru_cache(maxsize=1000000)
def _parse_prefix_cached(source_prefix: str, salt: int, profile_name: str):
    return parse_and_validate_prefix(source_prefix, salt, PROFILES[profile_name])


def nominal_schedule(profile: Profile) -> tuple[int, ...]:
    """Nominal source-rate schedule over the 68-nt main region only."""
    S = profile.source_bits
    n0 = len(ANCHOR_POSITIONS)
    usable = MAIN_NT - n0
    # n1 + n2 = usable; n1 + 2*n2 = S
    n_two = S - usable
    n_one = 2*usable - S
    if n_two < 0 or n_one < 0:
        raise ValueError("Profile source bits incompatible with main-region capacity")
    rates = [2] * MAIN_NT
    for p in ANCHOR_POSITIONS:
        rates[p] = 0
    ones = set()
    for a in ANCHOR_POSITIONS:
        for d in (-2, -1, 1, 2):
            p = a + d
            if 0 <= p < MAIN_NT and rates[p] == 2:
                ones.add(p)
                if len(ones) >= n_one:
                    break
        if len(ones) >= n_one:
            break
    while len(ones) < n_one:
        candidates = [i for i in range(MAIN_NT) if rates[i] == 2 and i not in ones]
        marked = list(ANCHOR_POSITIONS) + list(ones)
        p = max(candidates, key=lambda i: (min(abs(i-j) for j in marked), i))
        ones.add(p)
    for p in ones:
        rates[p] = 1
    assert rates.count(0) == n0 and rates.count(1) == n_one and rates.count(2) == n_two
    assert sum(rates) == S
    return tuple(rates)


@lru_cache(maxsize=32)
def _nominal_schedule_cached(profile_name: str) -> tuple[int, ...]:
    return nominal_schedule(PROFILES[profile_name])


def _gc_feasible_after(new_gc: int, position: int) -> bool:
    used = position + 1
    rem = BODY_NT - used
    lo = math.ceil(0.40 * BODY_NT)
    hi = math.floor(0.60 * BODY_NT)
    return (new_gc <= hi) and (new_gc + rem >= lo)


def allowed_bases(state: BioState, position: int) -> tuple[str, ...]:
    out = []
    for b in BASES:
        if b == state.last_base and state.run_len >= 3:
            continue
        new_gc = state.gc_count + int(b in "GC")
        if not _gc_feasible_after(new_gc, position):
            continue
        out.append(b)
    return tuple(out)


def _capacity_from_allowed(n: int) -> int:
    return 2 if n >= 4 else 1 if n >= 2 else 0 if n >= 1 else -1


def _rate_for_state(profile: Profile, position: int, source_len: int, debt: int, state: BioState):
    if position >= MAIN_NT:
        return -1, debt, tuple()
    if source_len >= profile.source_bits:
        return 0, debt, allowed_bases(state, position)
    nominal = _nominal_schedule_cached(profile.name)[position]
    allowed = allowed_bases(state, position)
    cap = _capacity_from_allowed(len(allowed))
    if cap < 0:
        return -1, debt, allowed
    if position in ANCHOR_POSITIONS:
        r = 0
    else:
        desired = nominal
        if debt > 0 and nominal < 2:
            desired = min(2, nominal + debt)
        r = min(desired, cap, profile.source_bits-source_len)
    new_debt = max(0, debt + nominal - r)
    rem_bits = profile.source_bits - (source_len+r)
    rem_positions = [p for p in range(position+1, MAIN_NT) if p not in ANCHOR_POSITIONS]
    if rem_bits > 2*len(rem_positions):
        return -1, new_debt, allowed
    return r, new_debt, allowed


def _tie_hash(seed_id: int, salt: int, position: int, history: str, dna_suffix: str, base: str) -> int:
    material = (
        TIE_DOMAIN + int(seed_id).to_bytes(2, "big") + int(salt).to_bytes(8, "big")
        + int(position).to_bytes(2, "big") + bits_to_bytes(history[-HISTORY_BITS:])
        + dna_suffix.encode() + base.encode()
    )
    return int.from_bytes(blake2s(material, digest_size=8).digest(), "big")


@lru_cache(maxsize=800000)
def choose_mapping(profile_name: str, seed_id: int, salt: int, position: int, source_len: int,
                   debt: int, history: str, dna_suffix: str, last_base: str, run_len: int, gc_count: int):
    profile = PROFILES[profile_name]
    state = BioState(dna_suffix, last_base, run_len, gc_count)
    r, new_debt, allowed = _rate_for_state(profile, position, source_len, debt, state)
    if r < 0 or not allowed:
        return -1, new_debt, tuple()
    need = 1 << r
    ranked = sorted(
        allowed,
        key=lambda b: _tie_hash(seed_id, salt, position, history, dna_suffix, b),
        reverse=True,
    )
    return r, new_debt, tuple(ranked[:need])



# =============================================================================
# V5.2 FILE-ENCODING RESCUE MODE
# =============================================================================
# The frozen V5.2 mapper remains unchanged for normal codewords (seed 0..15).
# If all 16 normal seeds reach a mapping dead-end, the application uses a
# deterministic rescue mapper with seed IDs 16..31.  Rescue mode keeps the
# same 103-bit payload, checkpoints, rate schedule, tie hash, and 76-nt BODY,
# but relaxes per-position GC/homopolymer restrictions so a mapping-complete
# codeword always exists.  The decoder mirrors this rule and searches both
# fixed hypothesis sets; no original payload or hidden selected seed is needed.

RESCUE_SEED_START = DEFAULT_SEED_COUNT
RESCUE_SEED_COUNT = 16
RAW_DECODE_HYPOTHESES = DEFAULT_SEED_COUNT + RESCUE_SEED_COUNT


def _rate_for_state_rescue(
    profile: Profile,
    position: int,
    source_len: int,
    debt: int,
    state: BioState,
):
    """Capacity-safe V5.2 rescue rate with unrestricted A/C/G/T choices."""
    if position >= MAIN_NT:
        return -1, debt, tuple()
    allowed = tuple(BASES)
    if source_len >= profile.source_bits:
        return 0, debt, allowed

    nominal = _nominal_schedule_cached(profile.name)[position]
    if position in ANCHOR_POSITIONS:
        r = 0
    else:
        desired = nominal
        if debt > 0 and nominal < 2:
            desired = min(2, nominal + debt)
        r = min(desired, 2, profile.source_bits - source_len)

    new_debt = max(0, debt + nominal - r)
    rem_bits = profile.source_bits - (source_len + r)
    rem_positions = [p for p in range(position + 1, MAIN_NT) if p not in ANCHOR_POSITIONS]
    if rem_bits > 2 * len(rem_positions):
        return -1, new_debt, allowed
    return r, new_debt, allowed


@lru_cache(maxsize=800000)
def choose_mapping_rescue(
    profile_name: str,
    seed_id: int,
    salt: int,
    position: int,
    source_len: int,
    debt: int,
    history: str,
    dna_suffix: str,
    last_base: str,
    run_len: int,
    gc_count: int,
):
    profile = PROFILES[profile_name]
    state = BioState(dna_suffix, last_base, run_len, gc_count)
    r, new_debt, allowed = _rate_for_state_rescue(
        profile, position, source_len, debt, state
    )
    if r < 0 or not allowed:
        return -1, new_debt, tuple()
    need = 1 << r
    ranked = sorted(
        allowed,
        key=lambda b: _tie_hash(
            seed_id, salt, position, history, dna_suffix, b
        ),
        reverse=True,
    )
    return r, new_debt, tuple(ranked[:need])


def _choose_mapping_mode(
    mapping_mode: str,
    profile_name: str,
    seed_id: int,
    salt: int,
    position: int,
    source_len: int,
    debt: int,
    history: str,
    dna_suffix: str,
    last_base: str,
    run_len: int,
    gc_count: int,
):
    if mapping_mode == "rescue":
        return choose_mapping_rescue(
            profile_name,
            seed_id,
            salt,
            position,
            source_len,
            debt,
            history,
            dna_suffix,
            last_base,
            run_len,
            gc_count,
        )
    return choose_mapping(
        profile_name,
        seed_id,
        salt,
        position,
        source_len,
        debt,
        history,
        dna_suffix,
        last_base,
        run_len,
        gc_count,
    )


@lru_cache(maxsize=800000)
def _terminal_rank(profile_name: str, seed_id: int, salt: int, tail_index: int,
                   source_bits: str, dna_suffix: str, last_base: str,
                   run_len: int, gc_count: int, base: str) -> int:
    material = (
        TERMINAL_DOMAIN + profile_name.encode()
        + int(seed_id).to_bytes(2, "big") + int(salt).to_bytes(8, "big")
        + int(tail_index).to_bytes(1, "big") + bits_to_bytes(source_bits)
        + dna_suffix.encode() + last_base.encode()
        + int(run_len).to_bytes(1, "big") + int(gc_count).to_bytes(1, "big")
        + base.encode()
    )
    return int.from_bytes(blake2s(material, digest_size=8).digest(), "big")


def terminal_base(profile_name: str, seed_id: int, salt: int, position: int,
                  source_bits: str, state: BioState) -> str:
    """Deterministic payload-dependent 0-bit terminal verifier symbol."""
    if position < MAIN_NT or position >= BODY_NT:
        raise ValueError("terminal_base called outside terminal tail")
    allowed = allowed_bases(state, position)
    if not allowed:
        raise RuntimeError("No biology-feasible terminal base")
    tail_index = position - MAIN_NT
    ranked = sorted(
        allowed,
        key=lambda b: _terminal_rank(
            profile_name, seed_id, salt, tail_index, source_bits,
            state.dna_suffix, state.last_base, state.run_len, state.gc_count, b,
        ),
        reverse=True,
    )
    return ranked[0]


def terminal_base_rescue(
    profile_name: str,
    seed_id: int,
    salt: int,
    position: int,
    source_bits: str,
    state: BioState,
) -> str:
    """Deterministic rescue terminal symbol with unrestricted A/C/G/T."""
    if position < MAIN_NT or position >= BODY_NT:
        raise ValueError("terminal_base_rescue called outside terminal tail")
    tail_index = position - MAIN_NT
    ranked = sorted(
        BASES,
        key=lambda b: _terminal_rank(
            profile_name,
            seed_id,
            salt,
            tail_index,
            source_bits,
            state.dna_suffix,
            state.last_base,
            state.run_len,
            state.gc_count,
            b,
        ),
        reverse=True,
    )
    return ranked[0]


def _terminal_base_mode(
    mapping_mode: str,
    profile_name: str,
    seed_id: int,
    salt: int,
    position: int,
    source_bits: str,
    state: BioState,
) -> str:
    if mapping_mode == "rescue":
        return terminal_base_rescue(
            profile_name, seed_id, salt, position, source_bits, state
        )
    return terminal_base(
        profile_name, seed_id, salt, position, source_bits, state
    )



@dataclass(frozen=True)
class EncodeResult:
    profile: str
    seed_id: int
    salt: int
    payload_bits: str
    source_bits: str
    dna: str
    body_bio: dict
    central_bio: dict
    rate_counts: tuple[int, int, int]


def encode_raw(payload: str, salt: int, si_seq: str, profile: Profile, seed_id: int) -> EncodeResult:
    source = build_source(payload, salt, profile)
    cursor = 0
    debt = 0
    history = ""
    st = initial_bio_from_si(si_seq)
    dna = []
    counts = [0,0,0]

    # 68-nt main HEDGES region: carries 103 payload + 4 internal check bits.
    for position in range(MAIN_NT):
        r, new_debt, mapping = choose_mapping(
            profile.name, seed_id, salt, position, cursor, debt, history,
            st.dna_suffix, st.last_base, st.run_len, st.gc_count,
        )
        if r < 0 or not mapping:
            raise RuntimeError(f"No feasible mapping at main position {position}")
        chunk = source[cursor:cursor+r]
        if len(chunk) != r:
            raise RuntimeError("Rate/source mismatch")
        value = int(chunk, 2) if r else 0
        b = mapping[value]
        dna.append(b)
        counts[r] += 1
        st = _update_bio(st, b)
        history = (history + chunk)[-HISTORY_BITS:]
        cursor += r
        debt = new_debt

    if cursor != profile.source_bits:
        raise RuntimeError(f"Could encode only {cursor}/{profile.source_bits} source bits before terminal tail")

    # 8-nt terminal verifier tail. These symbols carry no new payload bits.
    for position in range(MAIN_NT, BODY_NT):
        b = terminal_base(profile.name, seed_id, salt, position, source, st)
        dna.append(b)
        counts[0] += 1
        st = _update_bio(st, b)

    body = "".join(dna)
    bb = bio_stats(body)
    cb = bio_stats(si_seq + body)
    if not bb["overall_bio_pass"] or not cb["overall_bio_pass"]:
        raise RuntimeError("Hard biology failed")
    return EncodeResult(profile.name, seed_id, salt, payload, source, body, bb, cb, tuple(counts))




def encode_raw_rescue(
    payload: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_id: int,
) -> EncodeResult:
    """Guaranteed mapping-complete rescue codeword for one 103-bit payload.

    The source/checkpoint construction and nominal rate schedule are unchanged.
    Only the biological allowed-base restriction is relaxed.  The matching
    rescue decoder uses the exact same deterministic mapping.
    """
    source = build_source(payload, salt, profile)
    cursor = 0
    debt = 0
    history = ""
    st_bio = initial_bio_from_si(si_seq)
    dna = []
    counts = [0, 0, 0]

    for position in range(MAIN_NT):
        r, new_debt, mapping = choose_mapping_rescue(
            profile.name,
            seed_id,
            salt,
            position,
            cursor,
            debt,
            history,
            st_bio.dna_suffix,
            st_bio.last_base,
            st_bio.run_len,
            st_bio.gc_count,
        )
        if r < 0 or not mapping:
            raise RuntimeError(
                f"Rescue mapper unexpectedly failed at main position {position}"
            )
        chunk = source[cursor:cursor + r]
        if len(chunk) != r:
            raise RuntimeError("Rescue rate/source mismatch")
        value = int(chunk, 2) if r else 0
        base = mapping[value]
        dna.append(base)
        counts[r] += 1
        st_bio = _update_bio(st_bio, base)
        history = (history + chunk)[-HISTORY_BITS:]
        cursor += r
        debt = new_debt

    if cursor != profile.source_bits:
        raise RuntimeError(
            f"Rescue mapper encoded only {cursor}/{profile.source_bits} source bits"
        )

    for position in range(MAIN_NT, BODY_NT):
        base = terminal_base_rescue(
            profile.name, seed_id, salt, position, source, st_bio
        )
        dna.append(base)
        counts[0] += 1
        st_bio = _update_bio(st_bio, base)

    body = "".join(dna)
    bb = bio_stats(body)
    cb = bio_stats(si_seq + body)
    return EncodeResult(
        profile.name,
        seed_id,
        salt,
        payload,
        source,
        body,
        bb,
        cb,
        tuple(counts),
    )



def _bounded_levenshtein(a: str, b: str, max_cost: int) -> int:
    """Exact Levenshtein distance with a max-cost band and early exit.

    Returns max_cost+1 if the true distance is larger than max_cost.
    Intended for short deterministic suffixes (typically <= 8 nt).
    """
    if max_cost < 0:
        return max_cost + 1
    if abs(len(a) - len(b)) > max_cost:
        return max_cost + 1
    if not a:
        return len(b) if len(b) <= max_cost else max_cost + 1
    if not b:
        return len(a) if len(a) <= max_cost else max_cost + 1

    INF = max_cost + 1
    prev = [j if j <= max_cost else INF for j in range(len(b) + 1)]
    for i, ca in enumerate(a, 1):
        cur = [i if i <= max_cost else INF] + [INF] * len(b)
        j0 = max(1, i - max_cost)
        j1 = min(len(b), i + max_cost)
        row_min = INF
        for j in range(j0, j1 + 1):
            sub = prev[j-1] + (ca != b[j-1])
            delete = prev[j] + 1
            insert = cur[j-1] + 1
            v = min(sub, delete, insert)
            cur[j] = v
            if v < row_min:
                row_min = v
        if row_min > max_cost:
            return INF
        prev = cur
    d = prev[len(b)]
    return d if d <= max_cost else INF


def _deterministic_suffix(profile: Profile, seed_id: int, salt: int,
                          code_pos: int, source_bits: str, history: str,
                          state: BioState, debt: int,
                          mapping_mode: str = "strict") -> str | None:
    """Generate all remaining code symbols once source_bits is complete.

    Because no new information remains, every remaining main-region position
    has rate 0 and the terminal verifier is also deterministic. Returning None
    means the biological/mapping constraints make this hypothesis impossible.
    """
    if len(source_bits) != profile.source_bits:
        return None
    if not (0 <= code_pos <= BODY_NT):
        return None

    st = state
    hist = history
    d = debt
    out = []

    for position in range(code_pos, min(MAIN_NT, BODY_NT)):
        r, new_debt, mapping = _choose_mapping_mode(
            mapping_mode,
            profile.name, seed_id, salt, position, len(source_bits), d, hist,
            st.dna_suffix, st.last_base, st.run_len, st.gc_count,
        )
        if r != 0 or len(mapping) != 1:
            return None
        b = mapping[0]
        out.append(b)
        st = _update_bio(st, b)
        d = new_debt

    for position in range(max(code_pos, MAIN_NT), BODY_NT):
        try:
            b = _terminal_base_mode(
                mapping_mode, profile.name, seed_id, salt, position, source_bits, st
            )
        except Exception:
            return None
        out.append(b)
        st = _update_bio(st, b)

    return "".join(out)


@dataclass(frozen=True)
class BodyState:
    code_pos: int
    recv_pos: int
    source_bits: str
    history: str
    bio: BioState
    debt: int
    edits: int


@dataclass(frozen=True)
class RawCandidate:
    payload_bits: str
    seed_id: int
    edit_distance: int


@dataclass(frozen=True)
class DecodeResult:
    candidates: tuple[RawCandidate, ...]
    states_expanded: int
    seconds: float


def _length_lb(s: BodyState, received_len: int) -> int:
    return abs((BODY_NT-s.code_pos) - (received_len-s.recv_pos))


def body_list_decode_raw(received: str, salt: int, si_seq: str, profile: Profile,
                         seed_id: int, max_edits: int,
                         mapping_mode: str = "strict") -> DecodeResult:
    """Exact bounded IDS list decoder with deterministic-suffix closure.

    The search is exhaustive over every main-region hypothesis within
    max_edits. As soon as a hypothesis has reconstructed the complete source,
    the rest of its codeword is deterministic. We then compute the exact
    remaining edit distance in one bounded Levenshtein DP instead of branching
    through every S/I/D realization of the suffix.
    """
    if any(b not in BASES for b in received):
        raise ValueError("DNA must contain only A/C/G/T")
    t0 = time.perf_counter()
    R = len(received)
    start = BodyState(0,0,"","",initial_bio_from_si(si_seq),0,0)
    buckets = defaultdict(list)
    buckets[0].append(start)
    best_seen = {}
    completions = {}
    states = 0

    for cost in range(max_edits+1):
        frontier = buckets[cost]
        idx = 0
        while idx < len(frontier):
            s = frontier[idx]; idx += 1
            if s.edits > max_edits or s.code_pos > BODY_NT or s.recv_pos > R:
                continue
            if s.edits + _length_lb(s,R) > max_edits:
                continue
            key = (s.code_pos, s.recv_pos, s.source_bits)
            prev_best = best_seen.get(key)
            if prev_best is not None and prev_best <= s.edits:
                continue
            best_seen[key] = s.edits

            ok, payload_prefix = _parse_prefix_cached(s.source_bits, salt, profile.name)
            if not ok:
                continue

            # V5.1 exact closure: once the full protected source is known, all
            # remaining code symbols are deterministic, including any 0-bit
            # main positions and the complete terminal tail.
            if len(s.source_bits) == profile.source_bits:
                if len(payload_prefix) != profile.payload_bits:
                    continue
                if build_source(payload_prefix, salt, profile) != s.source_bits:
                    continue
                suffix = _deterministic_suffix(
                    profile, seed_id, salt, s.code_pos, s.source_bits,
                    s.history, s.bio, s.debt, mapping_mode,
                )
                if suffix is None:
                    continue
                rem_budget = max_edits - s.edits
                extra = _bounded_levenshtein(suffix, received[s.recv_pos:], rem_budget)
                if extra <= rem_budget:
                    final_edits = s.edits + extra
                    old = completions.get(payload_prefix)
                    if old is None or final_edits < old:
                        completions[payload_prefix] = final_edits
                continue

            if s.code_pos == BODY_NT:
                continue

            states += 1

            # Source is still incomplete, so we must be inside the main HEDGES
            # region and branch exactly as in V5.
            if s.code_pos >= MAIN_NT:
                continue

            r, new_debt, mapping = _choose_mapping_mode(
                mapping_mode,
                profile.name, seed_id, salt, s.code_pos, len(s.source_bits), s.debt, s.history,
                s.bio.dna_suffix, s.bio.last_base, s.bio.run_len, s.bio.gc_count,
            )
            if r < 0 or not mapping:
                continue

            # Received insertion: consume a read base, emit no code base.
            if s.recv_pos < R and s.edits < max_edits:
                ns = BodyState(
                    s.code_pos, s.recv_pos+1, s.source_bits, s.history,
                    s.bio, s.debt, s.edits+1
                )
                if ns.edits + _length_lb(ns,R) <= max_edits:
                    buckets[ns.edits].append(ns)

            for value, predicted in enumerate(mapping):
                chunk = f"{value:0{r}b}" if r else ""
                new_source = s.source_bits + chunk
                if len(new_source) > profile.source_bits:
                    continue
                ok2, _ = _parse_prefix_cached(new_source, salt, profile.name)
                if not ok2:
                    continue
                nh = (s.history + chunk)[-HISTORY_BITS:]
                nb = _update_bio(s.bio, predicted)

                # Match/substitution.
                if s.recv_pos < R:
                    nc = s.edits + int(predicted != received[s.recv_pos])
                    if nc <= max_edits:
                        ns = BodyState(
                            s.code_pos+1, s.recv_pos+1, new_source, nh, nb,
                            new_debt, nc
                        )
                        if ns.edits + _length_lb(ns,R) <= max_edits:
                            buckets[nc].append(ns)

                # Code deletion: emit predicted code symbol but consume no read base.
                if s.edits < max_edits:
                    ns = BodyState(
                        s.code_pos+1, s.recv_pos, new_source, nh, nb,
                        new_debt, s.edits+1
                    )
                    if ns.edits + _length_lb(ns,R) <= max_edits:
                        buckets[ns.edits].append(ns)

    cands = tuple(
        RawCandidate(p, seed_id, d)
        for p,d in sorted(completions.items(), key=lambda x:(x[1],x[0]))
    )
    return DecodeResult(cands, states, time.perf_counter()-t0)


def raw_multiseed_decode(received: str, salt: int, si_seq: str, profile: Profile,
                         seed_count: int, max_edits: int) -> DecodeResult:
    t0 = time.perf_counter()
    best = {}
    states = 0
    for seed_id in range(seed_count):
        r = body_list_decode_raw(received, salt, si_seq, profile, seed_id, max_edits)
        states += r.states_expanded
        for c in r.candidates:
            k = (c.payload_bits, c.seed_id)
            old = best.get(k)
            if old is None or c.edit_distance < old:
                best[k] = c.edit_distance
    cands = tuple(
        RawCandidate(p,s,d) for (p,s),d in sorted(best.items(), key=lambda x:(x[1],x[0][1],x[0][0]))
    )
    return DecodeResult(cands, states, time.perf_counter()-t0)



def raw_all_modes_decode(
    received: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    max_edits: int,
) -> DecodeResult:
    """Fail-closed raw list decode across 16 strict + 16 rescue hypotheses."""
    t0 = time.perf_counter()
    best = {}
    states = 0

    strict = raw_multiseed_decode(
        received,
        salt,
        si_seq,
        profile,
        DEFAULT_SEED_COUNT,
        max_edits,
    )
    states += strict.states_expanded
    for c in strict.candidates:
        key = (c.payload_bits, c.seed_id)
        old = best.get(key)
        if old is None or c.edit_distance < old:
            best[key] = c.edit_distance

    for seed_id in range(RESCUE_SEED_START, RAW_DECODE_HYPOTHESES):
        r = body_list_decode_raw(
            received,
            salt,
            si_seq,
            profile,
            seed_id,
            max_edits,
            mapping_mode="rescue",
        )
        states += r.states_expanded
        for c in r.candidates:
            key = (c.payload_bits, c.seed_id)
            old = best.get(key)
            if old is None or c.edit_distance < old:
                best[key] = c.edit_distance

    cands = tuple(
        RawCandidate(payload, seed_id, dist)
        for (payload, seed_id), dist in sorted(
            best.items(),
            key=lambda x: (x[1], x[0][1], x[0][0]),
        )
    )
    return DecodeResult(cands, states, time.perf_counter() - t0)


@dataclass(frozen=True)
class SeedCertificate:
    seed_id: int
    feasible: bool
    certified: bool
    competitor_payloads: int
    raw_candidate_pairs: int
    states: int
    seconds: float
    body_gc: float | None
    max_hp: int | None


@dataclass(frozen=True)
class SelectionResult:
    status: str
    encoding: EncodeResult | None
    certificates: tuple[SeedCertificate, ...]


def _selection_tie(payload: str, salt: int, profile_name: str, seed_id: int) -> int:
    material = (
        SELECT_DOMAIN + profile_name.encode() + int(salt).to_bytes(8,"big")
        + int(seed_id).to_bytes(2,"big") + bits_to_bytes(payload)
    )
    return int.from_bytes(blake2s(material, digest_size=8).digest(), "big")


@lru_cache(maxsize=50000)
def select_verified_encoding(payload: str, salt: int, si_seq: str, profile_name: str,
                             seed_count: int = DEFAULT_SEED_COUNT, cert_radius: int = DEFAULT_CERT_RADIUS,
                             strict: bool = True) -> SelectionResult:
    profile = PROFILES[profile_name]
    certs = []
    encodings = {}
    for seed_id in range(seed_count):
        try:
            enc = encode_raw(payload, salt, si_seq, profile, seed_id)
        except Exception:
            certs.append(SeedCertificate(seed_id, False, False, 10**9, 0, 0, 0.0, None, None))
            continue
        encodings[seed_id] = enc
        dec = raw_multiseed_decode(enc.dna, salt, si_seq, profile, seed_count, cert_radius)
        competitor_payloads = len({c.payload_bits for c in dec.candidates if c.payload_bits != payload})
        certs.append(SeedCertificate(
            seed_id=seed_id,
            feasible=True,
            certified=(competitor_payloads == 0),
            competitor_payloads=competitor_payloads,
            raw_candidate_pairs=len(dec.candidates),
            states=dec.states_expanded,
            seconds=dec.seconds,
            body_gc=enc.body_bio["gc_percent"],
            max_hp=enc.body_bio["longest_homopolymer"],
        ))

    feasible = [c for c in certs if c.feasible]
    certified = [c for c in feasible if c.certified]
    pool = certified if certified else ([] if strict else feasible)
    if not pool:
        return SelectionResult("ENCODE_UNCERTIFIED" if feasible else "ENCODE_FAIL", None, tuple(certs))

    def score(c: SeedCertificate):
        enc = encodings[c.seed_id]
        gc_slack = min(enc.body_bio["gc_percent"]-40.0, 60.0-enc.body_bio["gc_percent"])
        return (
            int(c.certified),
            -c.competitor_payloads,
            gc_slack,
            -c.max_hp,
            -c.states,
            _selection_tie(payload, salt, profile_name, c.seed_id),
        )

    chosen = max(pool, key=score)
    return SelectionResult("OK" if chosen.certified else "OK_UNCERTIFIED", encodings[chosen.seed_id], tuple(certs))


@dataclass(frozen=True)
class VerifiedCandidate:
    payload_bits: str
    seed_id: int
    edit_distance: int


@dataclass(frozen=True)
class VerifiedDecodeResult:
    status: str
    candidate: VerifiedCandidate | None
    candidates: tuple[VerifiedCandidate, ...]
    raw_candidate_pairs: int
    states_expanded: int
    seconds: float


def verified_decode(received: str, salt: int, si_seq: str, profile: Profile,
                    seed_count: int = DEFAULT_SEED_COUNT, max_edits: int = 1,
                    cert_radius: int = DEFAULT_CERT_RADIUS, strict: bool = True) -> VerifiedDecodeResult:
    t0 = time.perf_counter()
    raw = raw_multiseed_decode(received, salt, si_seq, profile, seed_count, max_edits)
    accepted = {}
    for c in raw.candidates:
        sel = select_verified_encoding(c.payload_bits, salt, si_seq, profile.name,
                                       seed_count, cert_radius, strict)
        if sel.encoding is None:
            continue
        if sel.encoding.seed_id != c.seed_id:
            continue
        old = accepted.get(c.payload_bits)
        vc = VerifiedCandidate(c.payload_bits, c.seed_id, c.edit_distance)
        if old is None or vc.edit_distance < old.edit_distance:
            accepted[c.payload_bits] = vc
    cands = tuple(sorted(accepted.values(), key=lambda c:(c.edit_distance,c.seed_id,c.payload_bits)))
    if len(cands) == 0:
        status, cand = "FAIL", None
    elif len(cands) == 1:
        status, cand = "OK", cands[0]
    else:
        status, cand = "AMBIGUOUS", None
    return VerifiedDecodeResult(status, cand, cands, len(raw.candidates), raw.states_expanded, time.perf_counter()-t0)


def fixed_seed_decode(received: str, salt: int, si_seq: str, profile: Profile,
                      seed_id: int = 0, max_edits: int = 1):
    raw = body_list_decode_raw(received, salt, si_seq, profile, seed_id, max_edits)
    by_payload = {}
    for c in raw.candidates:
        old = by_payload.get(c.payload_bits)
        if old is None or c.edit_distance < old.edit_distance:
            by_payload[c.payload_bits] = c
    cands = tuple(sorted(by_payload.values(), key=lambda c:(c.edit_distance,c.payload_bits)))
    if len(cands) == 0:
        status, cand = "FAIL", None
    elif len(cands) == 1:
        status, cand = "OK", cands[0]
    else:
        status, cand = "AMBIGUOUS", None
    return status, cand, cands, raw.states_expanded, raw.seconds


# ---------------- SI list decode / optional joint decoding ----------------
@dataclass(frozen=True)
class SICandidate:
    strand_id: int
    si_sequence: str
    split: int
    edit_distance: int


def si_list_candidates(central_read: str, codebook: Dict[int,str], max_si_edits: int):
    found = {}
    for split in range(max(1,SI_NT-max_si_edits), min(len(central_read),SI_NT+max_si_edits)+1):
        obs = central_read[:split]
        for strand_id, codeword in codebook.items():
            if abs(len(obs)-len(codeword)) > max_si_edits:
                continue
            d = levenshtein(obs, codeword)
            if d <= max_si_edits:
                found[(strand_id,split)] = SICandidate(strand_id,codeword,split,d)
    return tuple(sorted(found.values(), key=lambda x:(x.edit_distance,x.strand_id,x.split)))


@dataclass(frozen=True)
class JointCandidate:
    strand_id: int
    payload_bits: str
    seed_id: int
    si_edit_distance: int
    body_edit_distance: int


def joint_verified_decode(central_read: str, codebook: Dict[int,str], profile: Profile,
                          seed_count=DEFAULT_SEED_COUNT, max_si_edits=1, max_body_edits=1,
                          cert_radius=DEFAULT_CERT_RADIUS, strict=True):
    si_cands = si_list_candidates(central_read, codebook, max_si_edits)
    out = {}
    total_states = 0
    t0 = time.perf_counter()
    for sc in si_cands:
        body_read = central_read[sc.split:]
        if abs(len(body_read)-BODY_NT) > max_body_edits:
            continue
        r = verified_decode(body_read, sc.strand_id, sc.si_sequence, profile,
                            seed_count, max_body_edits, cert_radius, strict)
        total_states += r.states_expanded
        for bc in r.candidates:
            k = (sc.strand_id, bc.payload_bits)
            jc = JointCandidate(sc.strand_id,bc.payload_bits,bc.seed_id,sc.edit_distance,bc.edit_distance)
            old = out.get(k)
            if old is None or (jc.si_edit_distance+jc.body_edit_distance) < (old.si_edit_distance+old.body_edit_distance):
                out[k] = jc
    cands = tuple(sorted(out.values(), key=lambda c:(c.si_edit_distance+c.body_edit_distance,c.strand_id,c.payload_bits)))
    status = "OK" if len(cands)==1 else "FAIL" if len(cands)==0 else "AMBIGUOUS"
    return status, (cands[0] if len(cands)==1 else None), cands, len(si_cands), total_states, time.perf_counter()-t0




def architecture_summary(profile_name: str = "P103_T8_C4") -> dict:
    p = PROFILES[profile_name]
    sched = nominal_schedule(p)
    return {
        "profile": p.name,
        "payload_bits": p.payload_bits,
        "internal_check_bits": p.check_bits,
        "source_bits": p.source_bits,
        "body_nt": BODY_NT,
        "main_nt": MAIN_NT,
        "terminal_nt": TERMINAL_NT,
        "anchors": ANCHOR_POSITIONS,
        "nominal_0_main": sched.count(0),
        "nominal_1_main": sched.count(1),
        "nominal_2_main": sched.count(2),
        "body_density": p.body_density,
        "full125_density": p.full125_density,
    }

@dataclass(frozen=True)
class ErrorTrace:
    substitutions: tuple[int,...]
    deletions: tuple[int,...]
    insertion_slots: tuple[int,...]


def inject_exact_errors(clean: str, n_sub: int, n_del: int, n_ins: int, rng: random.Random):
    if n_sub+n_del > len(clean):
        raise ValueError("Too many sub+del")
    idx = list(range(len(clean)))
    rng.shuffle(idx)
    dels = set(idx[:n_del])
    subs = set(idx[n_del:n_del+n_sub])
    chars=[]
    for i,b in enumerate(clean):
        if i in dels:
            continue
        if i in subs:
            b = rng.choice([x for x in BASES if x != b])
        chars.append(b)
    slots=[]
    for _ in range(n_ins):
        slot=rng.randrange(len(chars)+1)
        chars.insert(slot,rng.choice(BASES)); slots.append(slot)
    return "".join(chars), ErrorTrace(tuple(sorted(subs)),tuple(sorted(dels)),tuple(sorted(slots)))

# =============================================================================
# V5.1 OPTIONAL FAST DETERMINISTIC SEED SELECTOR
# =============================================================================
# This mode preserves the V5 rate/tail architecture but changes which feasible
# mapper seed is selected. It avoids recursive K1 re-certification during
# decoding. Keep COMPAT_VERIFIED for exact V5 codeword compatibility.

@dataclass(frozen=True)
class FastSelectionResult:
    status: str
    encoding: EncodeResult | None
    feasible_seed_ids: tuple[int, ...]


@lru_cache(maxsize=100000)
def select_fast_encoding(payload: str, salt: int, si_seq: str, profile_name: str,
                         seed_count: int = DEFAULT_SEED_COUNT) -> FastSelectionResult:
    profile = PROFILES[profile_name]
    encodings = []
    for seed_id in range(seed_count):
        try:
            enc = encode_raw(payload, salt, si_seq, profile, seed_id)
        except Exception:
            continue
        encodings.append(enc)
    if not encodings:
        return FastSelectionResult("ENCODE_FAIL", None, tuple())

    def score(enc: EncodeResult):
        gc = enc.body_bio["gc_percent"]
        gc_slack = min(gc - 40.0, 60.0 - gc)
        return (
            gc_slack,
            -enc.body_bio["longest_homopolymer"],
            _selection_tie(payload, salt, profile_name, enc.seed_id),
        )

    chosen = max(encodings, key=score)
    return FastSelectionResult(
        "OK", chosen, tuple(sorted(enc.seed_id for enc in encodings))
    )


def fast_verified_decode(received: str, salt: int, si_seq: str, profile: Profile,
                         seed_count: int = DEFAULT_SEED_COUNT, max_edits: int = 1) -> VerifiedDecodeResult:
    """Fail-closed list decode using the cheap deterministic feasible-seed selector.

    Unlike COMPAT_VERIFIED, this does not recursively run K1 certification for
    each payload candidate. Each candidate payload is re-encoded across seeds,
    the unique deterministic feasible seed is selected, and candidates whose
    decoding seed does not equal that selected seed are rejected.
    """
    t0 = time.perf_counter()
    raw = raw_multiseed_decode(received, salt, si_seq, profile, seed_count, max_edits)
    accepted = {}
    for c in raw.candidates:
        sel = select_fast_encoding(c.payload_bits, salt, si_seq, profile.name, seed_count)
        if sel.encoding is None or sel.encoding.seed_id != c.seed_id:
            continue
        vc = VerifiedCandidate(c.payload_bits, c.seed_id, c.edit_distance)
        old = accepted.get(c.payload_bits)
        if old is None or vc.edit_distance < old.edit_distance:
            accepted[c.payload_bits] = vc
    cands = tuple(sorted(accepted.values(), key=lambda c:(c.edit_distance,c.seed_id,c.payload_bits)))
    if len(cands) == 0:
        status, cand = "FAIL", None
    elif len(cands) == 1:
        status, cand = "OK", cands[0]
    else:
        status, cand = "AMBIGUOUS", None
    return VerifiedDecodeResult(status, cand, cands, len(raw.candidates), raw.states_expanded, time.perf_counter()-t0)


def decode_v51(received: str, salt: int, si_seq: str, profile: Profile,
               seed_count: int = DEFAULT_SEED_COUNT, max_edits: int = 1,
               selector_mode: str = "COMPAT_VERIFIED",
               cert_radius: int = DEFAULT_CERT_RADIUS,
               strict: bool = True) -> VerifiedDecodeResult:
    mode = selector_mode.upper()
    if mode == "FAST_FEASIBLE":
        return fast_verified_decode(received, salt, si_seq, profile, seed_count, max_edits)
    if mode == "COMPAT_VERIFIED":
        return verified_decode(received, salt, si_seq, profile, seed_count, max_edits,
                               cert_radius, strict)
    raise ValueError("selector_mode must be FAST_FEASIBLE or COMPAT_VERIFIED")


# =============================================================================
# V5.2 EXACT PARALLEL SEED DECODER
# =============================================================================
#
# Scientific invariant:
#   Each mapper seed defines an independent exact bounded list-decode problem.
#   Parallel execution changes only the order in which these independent
#   searches are executed. Candidate aggregation is identical to V5.1.
#
# Default policy:
#   - K <= 2 : sequential V5.1 decoder (process overhead is not worthwhile)
#   - K >= 3 : process-level parallelism over seed hypotheses
#
# No beam search, nearest-neighbor shortcut, candidate truncation, or heuristic
# seed elimination is used.

DEFAULT_MAX_WORKERS = 0       # 0 = auto
DEFAULT_PARALLEL_MIN_EDITS = 3


def _resolve_workers(seed_count: int, max_workers: int | None) -> int:
    """Choose a conservative worker count.

    Auto mode is capped at 4 to avoid excessive Colab memory/process overhead.
    On the common 2-vCPU Colab runtime this resolves to 2.
    """
    if max_workers is None or int(max_workers) <= 0:
        cpu = os.cpu_count() or 2
        return max(1, min(int(seed_count), int(cpu), 4))
    return max(1, min(int(seed_count), int(max_workers)))


def _seed_decode_worker(args):
    """Top-level pickle-safe worker for one exact mapper-seed search."""
    received, salt, si_seq, profile_name, seed_id, max_edits = args
    profile = PROFILES[profile_name]
    r = body_list_decode_raw(
        received,
        int(salt),
        si_seq,
        profile,
        int(seed_id),
        int(max_edits),
    )
    # Return only simple/pickle-safe values.
    cand_rows = tuple(
        (c.payload_bits, int(c.seed_id), int(c.edit_distance))
        for c in r.candidates
    )
    return int(seed_id), cand_rows, int(r.states_expanded), float(r.seconds)


def _aggregate_seed_results(seed_results, elapsed: float) -> DecodeResult:
    """Deterministically aggregate per-seed exact candidate lists."""
    best = {}
    states = 0
    for seed_id, cand_rows, seed_states, _seed_seconds in sorted(
        seed_results, key=lambda x: x[0]
    ):
        states += int(seed_states)
        for payload_bits, candidate_seed, edit_distance in cand_rows:
            key = (payload_bits, int(candidate_seed))
            old = best.get(key)
            if old is None or int(edit_distance) < old:
                best[key] = int(edit_distance)

    cands = tuple(
        RawCandidate(payload, seed_id, dist)
        for (payload, seed_id), dist in sorted(
            best.items(),
            key=lambda x: (x[1], x[0][1], x[0][0]),
        )
    )
    return DecodeResult(cands, states, float(elapsed))


def raw_multiseed_decode_parallel(
    received: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_count: int,
    max_edits: int,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
) -> DecodeResult:
    """Exact V5.1 raw list decode, parallelized only across mapper seeds.

    The returned candidate set and summed state count must be identical to
    `raw_multiseed_decode()` for the same inputs. Only wall-clock execution
    differs.
    """
    workers = _resolve_workers(seed_count, max_workers)
    if workers <= 1 or seed_count <= 1:
        return raw_multiseed_decode(
            received, salt, si_seq, profile, seed_count, max_edits
        )

    t0 = time.perf_counter()
    tasks = [
        (
            received,
            int(salt),
            si_seq,
            profile.name,
            int(seed_id),
            int(max_edits),
        )
        for seed_id in range(seed_count)
    ]

    # `fork` is efficient on Linux/Colab. Fall back to default context where
    # fork is unavailable (e.g. Windows).
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()

    results = []
    executor = ProcessPoolExecutor(max_workers=workers, mp_context=ctx)
    try:
        futures = [executor.submit(_seed_decode_worker, task) for task in tasks]
        for fut in as_completed(futures):
            results.append(fut.result())
    except BaseException:
        # Important for notebook timeouts/interrupts: do not wait for all queued
        # futures before propagating the exception.
        for fut in futures:
            fut.cancel()
        executor.shutdown(wait=False, cancel_futures=True)
        raise
    else:
        executor.shutdown(wait=True)

    return _aggregate_seed_results(results, time.perf_counter() - t0)


def _raw_candidate_signature(result: DecodeResult):
    return tuple(
        sorted(
            (c.payload_bits, int(c.seed_id), int(c.edit_distance))
            for c in result.candidates
        )
    )


def assert_parallel_equivalence(
    received: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_count: int = DEFAULT_SEED_COUNT,
    max_edits: int = 3,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
) -> dict:
    """Run sequential and parallel raw decoders and assert exact equivalence."""
    seq = raw_multiseed_decode(
        received, salt, si_seq, profile, seed_count, max_edits
    )
    par = raw_multiseed_decode_parallel(
        received,
        salt,
        si_seq,
        profile,
        seed_count,
        max_edits,
        max_workers=max_workers,
    )
    same_candidates = _raw_candidate_signature(seq) == _raw_candidate_signature(par)
    same_states = int(seq.states_expanded) == int(par.states_expanded)
    if not same_candidates or not same_states:
        raise AssertionError(
            f"V5.2 regression mismatch: candidates={same_candidates}, "
            f"states={same_states}"
        )
    return {
        "same_candidates": same_candidates,
        "same_states": same_states,
        "sequential_s": float(seq.seconds),
        "parallel_s": float(par.seconds),
        "speedup": (
            float(seq.seconds / par.seconds) if par.seconds > 0 else float("inf")
        ),
        "states": int(seq.states_expanded),
        "candidate_pairs": len(seq.candidates),
        "workers": _resolve_workers(seed_count, max_workers),
    }


def _accept_verified_candidates(
    raw: DecodeResult,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_count: int,
    cert_radius: int,
    strict: bool,
) -> tuple[VerifiedCandidate, ...]:
    """Same post-list verification rule used by V5/V5.1."""
    accepted = {}
    for c in raw.candidates:
        sel = select_verified_encoding(
            c.payload_bits,
            salt,
            si_seq,
            profile.name,
            seed_count,
            cert_radius,
            strict,
        )
        if sel.encoding is None:
            continue
        if sel.encoding.seed_id != c.seed_id:
            continue
        vc = VerifiedCandidate(c.payload_bits, c.seed_id, c.edit_distance)
        old = accepted.get(c.payload_bits)
        if old is None or vc.edit_distance < old.edit_distance:
            accepted[c.payload_bits] = vc
    return tuple(
        sorted(
            accepted.values(),
            key=lambda c: (c.edit_distance, c.seed_id, c.payload_bits),
        )
    )


def verified_decode_v52(
    received: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_count: int = DEFAULT_SEED_COUNT,
    max_edits: int = 1,
    cert_radius: int = DEFAULT_CERT_RADIUS,
    strict: bool = True,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
    parallel_min_edits: int = DEFAULT_PARALLEL_MIN_EDITS,
) -> VerifiedDecodeResult:
    """V5.2 exact verified decoder.

    For low K it delegates to the original V5.1 sequential decoder.
    For K >= parallel_min_edits it parallelizes only the independent seed
    searches, then applies the unchanged verified seed-selection rule.
    """
    if int(max_edits) < int(parallel_min_edits):
        return verified_decode(
            received,
            salt,
            si_seq,
            profile,
            seed_count,
            max_edits,
            cert_radius,
            strict,
        )

    t0 = time.perf_counter()
    raw = raw_multiseed_decode_parallel(
        received,
        salt,
        si_seq,
        profile,
        seed_count,
        max_edits,
        max_workers=max_workers,
    )
    cands = _accept_verified_candidates(
        raw,
        salt,
        si_seq,
        profile,
        seed_count,
        cert_radius,
        strict,
    )

    if len(cands) == 0:
        status, cand = "FAIL", None
    elif len(cands) == 1:
        status, cand = "OK", cands[0]
    else:
        status, cand = "AMBIGUOUS", None

    return VerifiedDecodeResult(
        status,
        cand,
        cands,
        len(raw.candidates),
        raw.states_expanded,
        time.perf_counter() - t0,
    )


def fast_verified_decode_v52(
    received: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_count: int = DEFAULT_SEED_COUNT,
    max_edits: int = 1,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
    parallel_min_edits: int = DEFAULT_PARALLEL_MIN_EDITS,
) -> VerifiedDecodeResult:
    """Experimental FAST_FEASIBLE selector + V5.2 parallel raw search."""
    if int(max_edits) < int(parallel_min_edits):
        return fast_verified_decode(
            received, salt, si_seq, profile, seed_count, max_edits
        )

    t0 = time.perf_counter()
    raw = raw_multiseed_decode_parallel(
        received,
        salt,
        si_seq,
        profile,
        seed_count,
        max_edits,
        max_workers=max_workers,
    )
    accepted = {}
    for c in raw.candidates:
        sel = select_fast_encoding(
            c.payload_bits, salt, si_seq, profile.name, seed_count
        )
        if sel.encoding is None or sel.encoding.seed_id != c.seed_id:
            continue
        vc = VerifiedCandidate(c.payload_bits, c.seed_id, c.edit_distance)
        old = accepted.get(c.payload_bits)
        if old is None or vc.edit_distance < old.edit_distance:
            accepted[c.payload_bits] = vc

    cands = tuple(
        sorted(
            accepted.values(),
            key=lambda c: (c.edit_distance, c.seed_id, c.payload_bits),
        )
    )
    if len(cands) == 0:
        status, cand = "FAIL", None
    elif len(cands) == 1:
        status, cand = "OK", cands[0]
    else:
        status, cand = "AMBIGUOUS", None

    return VerifiedDecodeResult(
        status,
        cand,
        cands,
        len(raw.candidates),
        raw.states_expanded,
        time.perf_counter() - t0,
    )


def decode_v52(
    received: str,
    salt: int,
    si_seq: str,
    profile: Profile,
    seed_count: int = DEFAULT_SEED_COUNT,
    max_edits: int = 1,
    selector_mode: str = "COMPAT_VERIFIED",
    cert_radius: int = DEFAULT_CERT_RADIUS,
    strict: bool = True,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
    parallel_min_edits: int = DEFAULT_PARALLEL_MIN_EDITS,
) -> VerifiedDecodeResult:
    mode = selector_mode.upper()
    if mode == "COMPAT_VERIFIED":
        return verified_decode_v52(
            received,
            salt,
            si_seq,
            profile,
            seed_count,
            max_edits,
            cert_radius,
            strict,
            max_workers,
            parallel_min_edits,
        )
    if mode == "FAST_FEASIBLE":
        return fast_verified_decode_v52(
            received,
            salt,
            si_seq,
            profile,
            seed_count,
            max_edits,
            max_workers,
            parallel_min_edits,
        )
    raise ValueError(
        "selector_mode must be COMPAT_VERIFIED or FAST_FEASIBLE"
    )


def joint_verified_decode_v52(
    central_read: str,
    codebook: Dict[int, str],
    profile: Profile,
    seed_count: int = DEFAULT_SEED_COUNT,
    max_si_edits: int = 1,
    max_body_edits: int = 1,
    cert_radius: int = DEFAULT_CERT_RADIUS,
    strict: bool = True,
    max_workers: int | None = DEFAULT_MAX_WORKERS,
    parallel_min_edits: int = DEFAULT_PARALLEL_MIN_EDITS,
):
    """Joint noisy-SI + BODY decode using the V5.2 BODY decoder."""
    si_cands = si_list_candidates(central_read, codebook, max_si_edits)
    out = {}
    total_states = 0
    t0 = time.perf_counter()

    for sc in si_cands:
        body_read = central_read[sc.split:]
        if abs(len(body_read) - BODY_NT) > max_body_edits:
            continue

        r = verified_decode_v52(
            body_read,
            sc.strand_id,
            sc.si_sequence,
            profile,
            seed_count=seed_count,
            max_edits=max_body_edits,
            cert_radius=cert_radius,
            strict=strict,
            max_workers=max_workers,
            parallel_min_edits=parallel_min_edits,
        )
        total_states += r.states_expanded

        for bc in r.candidates:
            key = (sc.strand_id, bc.payload_bits)
            jc = JointCandidate(
                sc.strand_id,
                bc.payload_bits,
                bc.seed_id,
                sc.edit_distance,
                bc.edit_distance,
            )
            old = out.get(key)
            if old is None or (
                jc.si_edit_distance + jc.body_edit_distance
            ) < (
                old.si_edit_distance + old.body_edit_distance
            ):
                out[key] = jc

    cands = tuple(
        sorted(
            out.values(),
            key=lambda c: (
                c.si_edit_distance + c.body_edit_distance,
                c.strand_id,
                c.payload_bits,
            ),
        )
    )
    status = "OK" if len(cands) == 1 else "FAIL" if len(cands) == 0 else "AMBIGUOUS"
    return (
        status,
        cands[0] if len(cands) == 1 else None,
        cands,
        len(si_cands),
        total_states,
        time.perf_counter() - t0,
    )
