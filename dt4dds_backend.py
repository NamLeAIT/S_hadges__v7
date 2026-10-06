"""Official DT4DDS molecular channel adapter.

The decoder never imports this module and never receives the evaluation SID map.
DT4DDS models molecule populations; ART is responsible for FASTQ bases/qualities.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from functools import lru_cache
from importlib import metadata
from pathlib import Path
import hashlib
import random
import sys


PAPER = "https://doi.org/10.1038/s41467-023-41729-1"
REPOSITORY = "https://github.com/fml-ethz/dt4dds"
PROFILE_NAME = "DT4DDS 1.1.0 · Twist synthesis + storage + Q5 PCR"


@dataclass(frozen=True)
class MolecularSettings:
    seed: int = 12345
    synthesis_coverage: int = 100
    synthesis_lognormal_sigma: float = 0.30
    storage_half_lives: float = 0.0
    pcr_cycles: int = 20
    pcr_efficiency: float = 0.95
    pcr_template_volume_ul: float = 1.0
    pcr_reaction_volume_ul: float = 20.0
    aliquot_molecules_per_strand: int = 100


def _imports():
    try:
        import dt4dds
        from rapidfuzz.distance import Levenshtein
    except ImportError as exc:
        if exc.name == "dt4dds":
            raise RuntimeError(
                "DT4DDS is not installed in the Python interpreter running Streamlit "
                f"({sys.executable}). Run install_windows.bat from this project, then "
                "start with run_windows.bat or: python -m streamlit run app.py."
            ) from exc
        raise RuntimeError(
            "DT4DDS molecular dependencies are incomplete. "
            f"Import failed for {exc.name or str(exc)}: {exc}. "
            "Run install_windows.bat on Windows; no simplified fallback was used."
        ) from exc
    version = getattr(dt4dds, "__VERSION__", None)
    if version is None:
        version = metadata.version("dt4dds")
    if version != "1.1.0":
        raise RuntimeError(f"This build validates dt4dds==1.1.0; found {version}.")
    return dt4dds, Levenshtein


def available():
    return dependency_status()[0]


def dependency_status():
    try:
        dt4dds, _ = _imports()
        return True, f"DT4DDS {dt4dds.__VERSION__} molecular backend ready"
    except (ImportError, RuntimeError, metadata.PackageNotFoundError) as exc:
        return False, str(exc)


def validate(rows, cfg: MolecularSettings):
    if not rows:
        raise ValueError("No designed strands were provided.")
    if not 1 <= cfg.synthesis_coverage <= 10000:
        raise ValueError("Synthesis sampling coverage must be 1..10,000 molecules/strand.")
    if not 0 <= cfg.synthesis_lognormal_sigma <= 3:
        raise ValueError("Synthesis lognormal sigma must be 0..3.")
    if not 0 <= cfg.storage_half_lives <= 10:
        raise ValueError("Storage exposure must be 0..10 half-lives.")
    if not 0 <= cfg.pcr_cycles <= 40 or not 0 <= cfg.pcr_efficiency <= 1:
        raise ValueError("Invalid PCR cycle count or efficiency.")
    if not 1 <= cfg.aliquot_molecules_per_strand <= 10000:
        raise ValueError("Physical aliquot must be 1..10,000 molecules/strand.")
    primers = {(str(r["FBR"]), str(r["RBR"])) for r in rows}
    if len(primers) != 1:
        raise ValueError("The DT4DDS workflow requires one shared primer pair per run.")
    for row in rows:
        if set(str(row["Full strand"]).upper()) - set("ACGT"):
            raise ValueError("DT4DDS input strands must contain only A/C/G/T.")


def _central(row):
    return str(row["SI"]) + str(row["Payload"]) + str(row.get("Filler", ""))


def _pool_records(pool):
    return [(str(seq), int(count)) for seq, count in pool if int(count) > 0]


def _editop_counts(editops):
    counts = Counter()
    for op in editops:
        if op.tag == "replace":
            counts["S"] += 1
        elif op.tag == "insert":
            counts["I"] += 1
        elif op.tag == "delete":
            counts["D"] += 1
    return counts


def _mapper(reference_by_sid, levenshtein):
    references = tuple(sorted((int(sid), str(seq)) for sid, seq in reference_by_sid.items()))

    @lru_cache(maxsize=50000)
    def assign(observed):
        distances = []
        for sid, reference in references:
            d = levenshtein.distance(reference, observed)
            distances.append((d, sid, reference))
        best_distance = min(x[0] for x in distances)
        tied = [x for x in distances if x[0] == best_distance]
        _, sid, reference = tied[0]
        errors = _editop_counts(levenshtein.editops(reference, observed))
        return sid, reference, errors, len(tied) > 1

    return assign


def _snapshot(stage, pool, reference_by_sid, levenshtein):
    assign = _mapper(reference_by_sid, levenshtein)
    totals = Counter()
    by_sid = Counter()
    ambiguous = 0
    examples = []
    for sequence, count in _pool_records(pool):
        sid, reference, errors, is_ambiguous = assign(sequence)
        by_sid[sid] += count
        ambiguous += count if is_ambiguous else 0
        totals["molecules"] += count
        totals["variants"] += 1
        totals["input_nt"] += len(reference) * count
        totals["sub"] += errors["S"] * count
        totals["ins"] += errors["I"] * count
        totals["del"] += errors["D"] * count
        if sequence == reference:
            totals["exact_molecules"] += count
        if len(examples) < 20 and sequence != reference:
            examples.append({"Stage": stage, "Strand ID": sid, "Before": reference,
                             "After": sequence, "Copies": count, "S": errors["S"],
                             "I": errors["I"], "D": errors["D"],
                             "Assignment ambiguous": is_ambiguous})
    errors = totals["sub"] + totals["ins"] + totals["del"]
    summary = {
        "Stage": stage,
        "sub": totals["sub"], "ins": totals["ins"], "del": totals["del"],
        "input_nt": totals["input_nt"], "molecules": totals["molecules"],
        "variants": totals["variants"], "ambiguous_molecules": ambiguous,
        "exact_molecules": totals["exact_molecules"],
        "observed_edit_rate": errors / totals["input_nt"] if totals["input_nt"] else None,
        "measurement": "abundance-weighted observed edit burden vs nearest designed strand",
    }
    return summary, by_sid, examples, assign


def _source_hashes(dt4dds):
    root = Path(dt4dds.__file__).resolve().parent
    wanted = [
        root / "settings" / "array_synthesis_twist.yaml",
        root / "settings" / "aging.yaml",
        root / "settings" / "pcr_Q5.yaml",
        root / "settings" / "properties.yaml",
    ]
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in wanted}


def manifest(cfg: MolecularSettings, dt4dds=None):
    if dt4dds is None:
        dt4dds, _ = _imports()
    return {
        "profile": PROFILE_NAME,
        "molecular_engine": "dt4dds",
        "dt4dds_version": getattr(dt4dds, "__VERSION__", "1.1.0"),
        "paper": PAPER,
        "repository": REPOSITORY,
        "configuration": asdict(cfg),
        "source_sha256": _source_hashes(dt4dds),
        "workflow": ["Twist array synthesis", "primer attachment", "aging/storage",
                     "Q5 PCR", "physical aliquot", "ART Illumina HS25"],
        "measurement_note": (
            "DT4DDS does not expose per-molecule mutation lineage. Stage S/I/D values are "
            "abundance-weighted minimum-edit burden versus the nearest designed strand, "
            "not counts of newly generated events. SID assignment is evaluation-only."
        ),
        "calibration_scope": (
            "Published DT4DDS process profiles plus explicit workflow assumptions; "
            "not yet calibrated to this laboratory or synthesis vendor."
        ),
    }


def run(rows, cfg: MolecularSettings, progress=None):
    """Run synthesis -> storage -> PCR -> physical aliquot.

    Returns a checkpoint consumed by :func:`channel_runs.finish`. The hidden
    ``_templates`` field contains evaluation-assigned molecular variants only;
    ART FASTQ headers remain neutral and the public decoder never receives it.
    """
    validate(rows, cfg)
    dt4dds, levenshtein = _imports()
    random.seed(cfg.seed)
    dt4dds.config.set_random_seed(cfg.seed)
    dt4dds.properties.AmplificationEfficiency._efficiency_cache = {}
    dt4dds.properties.set_property_settings(
        dt4dds.settings.defaults.SequenceProperties(
            efficiency_distribution="normal",
            efficiency_params={"loc": 1.0, "scale": 0.0051},
        )
    )
    if progress:
        progress("Synthesis", 0.05)

    central_refs = {int(r["Strand ID"]): _central(r) for r in rows}
    full_refs = {int(r["Strand ID"]): str(r["Full strand"]) for r in rows}
    fbr, rbr = next(iter({(str(r["FBR"]), str(r["RBR"])) for r in rows}))
    reverse_primer = str(dt4dds.datastructures.Seq(rbr).reverse_complement())

    synthesis_settings = dt4dds.settings.defaults.ArraySynthesis_Twist(
        oligo_distribution_type="lognormal",
        oligo_distribution_params={"mean": 0, "sigma": cfg.synthesis_lognormal_sigma},
    )
    synthesis = dt4dds.processes.ArraySynthesis(synthesis_settings)
    synthesis.process(list(central_refs.values()))
    pool = synthesis.sample_by_counts(cfg.synthesis_coverage * len(rows), volume=1)
    pool = dt4dds.generators.attach_primers_to_pool(pool, fbr, reverse_primer)
    pool.volume = 1
    # Report all molecular stages against the same full-strand denominator.
    # The array synthesis errors are in the central oligo; primers are attached
    # cleanly as the explicit workflow step used by the reference scenario.
    syn_summary, _, syn_examples, _ = _snapshot("Synthesis", pool, full_refs, levenshtein)
    if progress:
        progress("Storage", 0.30)
    if cfg.storage_half_lives > 0:
        aging = dt4dds.processes.Aging(
            dt4dds.settings.defaults.Aging(), n_halflives=cfg.storage_half_lives
        )
        pool = aging.process(pool)
        pool.volume = 1
    storage_summary, _, storage_examples, _ = _snapshot("Storage", pool, full_refs, levenshtein)

    if progress:
        progress("PCR", 0.50)
    pcr_settings = dt4dds.settings.defaults.PCR_Q5(
        primers=(fbr, reverse_primer),
        template_volume=cfg.pcr_template_volume_ul,
        volume=cfg.pcr_reaction_volume_ul,
        efficiency_mean=cfg.pcr_efficiency,
        n_cycles=cfg.pcr_cycles,
    )
    pool = dt4dds.processes.PCR(pcr_settings).process(pool)
    pcr_summary, pcr_by_sid, pcr_examples, _ = _snapshot("PCR", pool, full_refs, levenshtein)

    if progress:
        progress("Physical sampling", 0.72)
    requested = cfg.aliquot_molecules_per_strand * len(rows)
    requested = min(requested, int(pool.n_oligos))
    if requested < 1:
        raise RuntimeError("No molecules remain for physical sampling.")
    aliquot = pool.sample_by_counts(requested, accurate=True)
    aliquot.volume = 1
    sampling_summary, aliquot_by_sid, sampling_examples, assign = _snapshot(
        "Sampling", aliquot, full_refs, levenshtein
    )
    templates = []
    for sequence, copies in _pool_records(aliquot):
        sid, _, _, ambiguous = assign(sequence)
        templates.append({"sid": sid, "template": sequence, "copies": copies,
                          "evaluation_assignment_ambiguous": ambiguous})
    abundance = []
    for sid in sorted(full_refs):
        abundance.append({"Strand ID": sid,
                          "PCR molecules": int(pcr_by_sid[sid]),
                          "Aliquot molecules": int(aliquot_by_sid[sid]),
                          "Read units": 0})
    if progress:
        progress("Molecular pool ready", 0.80)
    return {
        "engine": "DT4DDS",
        "completed_stage": "Sampling",
        "settings": {"scope": "Full strand", "mode": "DT4DDS molecular channel + ART",
                     "seed": cfg.seed, "molecular": asdict(cfg)},
        "stages": [syn_summary, storage_summary, pcr_summary, sampling_summary],
        "stage_quality": [],
        "abundance": abundance,
        "audit": syn_examples + storage_examples + pcr_examples + sampling_examples,
        "trace": [],
        "preset": manifest(cfg, dt4dds),
        "_templates": templates,
    }
