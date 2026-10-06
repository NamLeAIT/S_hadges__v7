from __future__ import annotations

"""Structural-randomness metrics used by the accepted App1 analysis UI.

The explicit paper-supported models retained here are:
- basic sequence quality and normalized base entropy;
- translational active-particle MSD: MSD(n)=4 D n + V^2 n^2;
- inverse-Ising polarization gamma plus a transparent pseudo-likelihood estimate
  of lambda/J and h (the exact missing SI-S4 closed form is not claimed);
- 3-input/1-output alpha/beta balance.
"""

from dataclasses import asdict, dataclass
from math import sqrt
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

DNA_ALPHABET = frozenset("ACGT")
DNA_TO_BITS: Mapping[str, str] = {"A": "00", "C": "01", "G": "10", "T": "11"}
TRANSLATIONAL_VECTORS: Mapping[str, Tuple[float, float]] = {
    "A": (1.0, 0.0), "C": (0.0, 1.0), "G": (-1.0, 0.0), "T": (0.0, -1.0)
}
LOGIC_CONTEXTS: Tuple[str, ...] = tuple(f"{i:03b}" for i in range(8))


class DNASequenceError(ValueError):
    pass


def clean_dna(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    seq = "".join(sequence.upper().split())
    if not seq:
        raise DNASequenceError("DNA sequence is empty")
    invalid = sorted(set(seq) - DNA_ALPHABET)
    if invalid:
        raise DNASequenceError("DNA contains unsupported symbols: " + ", ".join(invalid))
    return seq


def dna_to_bits(sequence: str) -> str:
    return "".join(DNA_TO_BITS[b] for b in clean_dna(sequence))


def _longest_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    longest = current = 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def _homopolymer_segment_count(seq: str, minimum_length: int = 2) -> int:
    if not seq:
        return 0
    count = 0
    run = 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            run += 1
        else:
            count += int(run >= minimum_length)
            run = 1
    count += int(run >= minimum_length)
    return count


def _normalized_entropy_from_counts(counts: Sequence[int]) -> float:
    v = np.asarray(counts, dtype=float)
    total = float(v.sum())
    if total <= 0:
        return float("nan")
    p = v[v > 0] / total
    h = -float(np.sum(p * np.log2(p)))
    return h / np.log2(len(v))


def sequence_quality_metrics(sequence: str, *, original_bit_count: Optional[int] = None) -> Dict[str, float]:
    seq = clean_dna(sequence)
    n = len(seq)
    counts = {b: seq.count(b) for b in "ACGT"}
    gc = (counts["G"] + counts["C"]) / n
    result: Dict[str, float] = {
        "length_nt": float(n),
        "a_fraction": counts["A"] / n,
        "c_fraction": counts["C"] / n,
        "g_fraction": counts["G"] / n,
        "t_fraction": counts["T"] / n,
        "gc_fraction": gc,
        "gc_deviation_from_0_5": abs(gc - 0.5),
        "longest_homopolymer": float(_longest_homopolymer(seq)),
        "homopolymer_segments_ge_2": float(_homopolymer_segment_count(seq, 2)),
        "base_entropy_normalized": _normalized_entropy_from_counts([counts[b] for b in "ACGT"]),
    }
    if original_bit_count is not None:
        result["information_density_bits_per_nt"] = float(original_bit_count) / n
    return result


def _fit_nonnegative_two_term_model(n_values: np.ndarray, y_values: np.ndarray, *, linear_multiplier: float):
    n_values = np.asarray(n_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    design = np.column_stack([linear_multiplier * n_values, n_values ** 2])
    candidates = []
    coef, *_ = np.linalg.lstsq(design, y_values, rcond=None)
    if np.all(coef >= 0):
        candidates.append(coef)
    q = design[:, 1]
    candidates.append(np.array([0.0, max(0.0, float(np.dot(q, y_values) / max(np.dot(q, q), 1e-30)))]))
    lin = design[:, 0]
    candidates.append(np.array([max(0.0, float(np.dot(lin, y_values) / max(np.dot(lin, lin), 1e-30))), 0.0]))
    candidates.append(np.array([0.0, 0.0]))
    best = min(candidates, key=lambda c: float(np.sum((y_values - design @ c) ** 2)))
    fitted = design @ best
    rss = float(np.sum((y_values - fitted) ** 2))
    tss = float(np.sum((y_values - np.mean(y_values)) ** 2))
    r2 = 1.0 - rss / tss if tss > 0 else 1.0
    return float(best[0]), sqrt(max(0.0, float(best[1]))), r2, fitted


def translational_active_particle_metrics(sequence: str, *, max_lag: int = 20) -> Dict[str, Any]:
    seq = clean_dna(sequence)
    max_lag = int(min(max_lag, len(seq) - 1))
    if max_lag < 2:
        raise DNASequenceError("Sequence is too short for translational MSD")
    steps = np.asarray([TRANSLATIONAL_VECTORS[b] for b in seq], dtype=float)
    positions = np.vstack([np.zeros((1, 2)), np.cumsum(steps, axis=0)])
    lags = np.arange(1, max_lag + 1, dtype=int)
    starts = np.arange(0, len(seq) - max_lag + 1)
    msd = np.empty(max_lag, dtype=float)
    for idx, lag in enumerate(lags):
        displacement = positions[starts + lag] - positions[starts]
        msd[idx] = float(np.mean(np.sum(displacement ** 2, axis=1)))
    diffusion, velocity, r2, fitted = _fit_nonnegative_two_term_model(lags.astype(float), msd, linear_multiplier=4.0)
    return {
        "velocity_V_L_per_nt": velocity,
        "diffusion_D_L2_per_nt": diffusion,
        "fit_r_squared": r2,
        "max_lag_nt": max_lag,
        "lags_nt": lags.tolist(),
        "msd": msd.tolist(),
        "msd_fitted": fitted.tolist(),
    }


def bits_to_matrix(bits: str, *, shape: Optional[Tuple[int, int]] = None) -> np.ndarray:
    if not bits or any(x not in "01" for x in bits):
        raise ValueError("bits must contain only 0/1")
    values = np.fromiter((int(x) for x in bits), dtype=np.int8)
    if shape is not None:
        rows, cols = shape
        required = rows * cols
        if rows < 2 or cols < 3 or len(values) < required:
            raise ValueError("Supplied matrix shape is not compatible with the sequence")
        return values[:required].reshape(rows, cols)
    rows = int(np.floor(np.sqrt(len(values))))
    while rows >= 2:
        cols = len(values) // rows
        if cols >= 3:
            return values[: rows * cols].reshape(rows, cols)
        rows -= 1
    raise DNASequenceError("Not enough bits for a 2x3 matrix")


def _fit_inverse_ising_pseudolikelihood(spins: np.ndarray, *, max_iter: int = 50, tolerance: float = 1e-9, ridge: float = 1e-8):
    neighbour_sum = np.zeros_like(spins, dtype=float)
    degree = np.zeros_like(spins, dtype=float)
    neighbour_sum[1:, :] += spins[:-1, :]; degree[1:, :] += 1
    neighbour_sum[:-1, :] += spins[1:, :]; degree[:-1, :] += 1
    neighbour_sum[:, 1:] += spins[:, :-1]; degree[:, 1:] += 1
    neighbour_sum[:, :-1] += spins[:, 1:]; degree[:, :-1] += 1
    valid = degree > 0
    x = neighbour_sum[valid].astype(float)
    y = ((spins[valid] + 1) / 2).astype(float)
    design = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2, dtype=float)
    converged = False
    iteration = 0
    for iteration in range(1, max_iter + 1):
        linear = np.clip(design @ beta, -40.0, 40.0)
        p = 1.0 / (1.0 + np.exp(-linear))
        w = np.maximum(p * (1.0 - p), 1e-9)
        gradient = design.T @ (y - p) - ridge * beta
        hessian = design.T @ (w[:, None] * design) + ridge * np.eye(2)
        try:
            step = np.linalg.solve(hessian, gradient)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(hessian) @ gradient
        new = beta + step
        if np.linalg.norm(new - beta) < tolerance:
            beta = new; converged = True; break
        beta = new
    return float(beta[1] / 2), float(beta[0] / 2), converged, iteration


def inverse_ising_metrics(sequence: str, *, matrix_shape: Optional[Tuple[int, int]] = None) -> Dict[str, Any]:
    matrix = bits_to_matrix(dna_to_bits(sequence), shape=matrix_shape)
    spins = 2.0 * matrix.astype(float) - 1.0
    gamma = float(np.mean(spins))
    horizontal = spins[:, :-1] * spins[:, 1:]
    vertical = spins[:-1, :] * spins[1:, :]
    corr = float((horizontal.sum() + vertical.sum()) / (horizontal.size + vertical.size))
    interaction, field, converged, iterations = _fit_inverse_ising_pseudolikelihood(spins)
    return {
        "polarization_gamma": gamma,
        "abs_polarization_gamma": abs(gamma),
        "ising_lambda": interaction,
        "abs_ising_lambda": abs(interaction),
        "ising_h": field,
        "abs_ising_h": abs(field),
        "nearest_neighbour_correlation": corr,
        "optimizer_converged": converged,
        "optimizer_iterations": iterations,
        "estimation_method": "standard_pseudolikelihood_not_missing_paper_S4_formula",
    }


def logic_3input_1output_metrics(sequence: str, *, matrix_shape: Optional[Tuple[int, int]] = None, layout: str = "paper_grid") -> Dict[str, Any]:
    bits = dna_to_bits(sequence)
    zero = {c: 0 for c in LOGIC_CONTEXTS}; one = {c: 0 for c in LOGIC_CONTEXTS}
    if layout == "paper_grid":
        matrix = bits_to_matrix(bits, shape=matrix_shape)
        rows, cols = matrix.shape
        for r in range(1, rows):
            for c in range(2, cols):
                ctx = f"{matrix[r,c-2]}{matrix[r,c-1]}{matrix[r-1,c]}"
                (zero if int(matrix[r,c]) == 0 else one)[ctx] += 1
    elif layout == "linear":
        for i in range(len(bits) - 3):
            ctx = bits[i:i+3]
            (zero if bits[i+3] == "0" else one)[ctx] += 1
    else:
        raise ValueError("layout must be paper_grid or linear")
    totals = {c: zero[c] + one[c] for c in LOGIC_CONTEXTS}
    n = sum(totals.values())
    if n <= 0:
        raise DNASequenceError("Not enough data for 3-input/1-output analysis")
    alpha = {c: (zero[c] / totals[c] if totals[c] else None) for c in LOGIC_CONTEXTS}
    beta = {c: totals[c] / n for c in LOGIC_CONTEXTS}
    observed_alpha = np.asarray([v for v in alpha.values() if v is not None], dtype=float)
    all_beta = np.asarray(list(beta.values()), dtype=float)
    return {
        "alpha": alpha,
        "beta": beta,
        "alpha_rmse_from_0_5": float(np.sqrt(np.mean((observed_alpha - 0.5) ** 2))),
        "beta_rmse_from_0_125": float(np.sqrt(np.mean((all_beta - 0.125) ** 2))),
        "missing_alpha_contexts": [c for c, v in alpha.items() if v is None],
        "counts_zero": zero,
        "counts_one": one,
        "sample_count": int(n),
        "layout": layout,
    }


@dataclass
class DNASequenceRandomnessReport:
    sequence_quality: Dict[str, Any]
    translational_active_particle: Dict[str, Any]
    inverse_ising: Dict[str, Any]
    logic_3input_1output: Dict[str, Any]
    metadata: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def app_summary(self) -> Dict[str, Any]:
        q = self.sequence_quality; t = self.translational_active_particle; i = self.inverse_ising; l = self.logic_3input_1output
        out = {
            "DNA length (nt)": int(q["length_nt"]),
            "Base entropy normalized": q["base_entropy_normalized"],
            "GC fraction": q["gc_fraction"],
            "GC deviation from 0.5": q["gc_deviation_from_0_5"],
            "Longest homopolymer": int(q["longest_homopolymer"]),
            "Homopolymer segments (2+ bases)": int(q["homopolymer_segments_ge_2"]),
            "Translational velocity V": t["velocity_V_L_per_nt"],
            "Translational diffusion D": t["diffusion_D_L2_per_nt"],
            "|Polarization gamma|": i["abs_polarization_gamma"],
            "|Inverse-Ising lambda|": i["abs_ising_lambda"],
            "|Inverse-Ising h|": i["abs_ising_h"],
            "Alpha RMSE from 0.5": l["alpha_rmse_from_0_5"],
            "Beta RMSE from 0.125": l["beta_rmse_from_0_125"],
        }
        if "information_density_bits_per_nt" in q:
            out["Information density (bits/nt)"] = q["information_density_bits_per_nt"]
        return out


def analyze_dna_randomness(sequence: str, *, original_bit_count: Optional[int] = None, matrix_shape: Optional[Tuple[int, int]] = None, max_lag: int = 20, logic_layout: str = "paper_grid") -> DNASequenceRandomnessReport:
    seq = clean_dna(sequence)
    return DNASequenceRandomnessReport(
        sequence_quality=sequence_quality_metrics(seq, original_bit_count=original_bit_count),
        translational_active_particle=translational_active_particle_metrics(seq, max_lag=max_lag),
        inverse_ising=inverse_ising_metrics(seq, matrix_shape=matrix_shape),
        logic_3input_1output=logic_3input_1output_metrics(seq, matrix_shape=matrix_shape, layout=logic_layout),
        metadata={
            "logic_layout": logic_layout,
            "max_lag": min(max_lag, len(seq)-1),
            "inverse_ising_note": "lambda/h are pseudo-likelihood estimates because the uploaded manuscript's exact SI-S4 expression was unavailable",
        },
    )


def compare_sequences(sequences: Mapping[str, str], *, original_bit_counts: Optional[Mapping[str, int]] = None, max_lag: int = 20) -> List[Dict[str, Any]]:
    rows = []
    for name, seq in sequences.items():
        bit_count = original_bit_counts.get(name) if original_bit_counts else None
        rep = analyze_dna_randomness(seq, original_bit_count=bit_count, max_lag=max_lag)
        row = {"Sequence": name}; row.update(rep.app_summary()); rows.append(row)
    return rows
