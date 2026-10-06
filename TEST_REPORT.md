# v13.1 compression quality dashboard — 2026-10-06

Targeted regression suite: **11 passed, 1 skipped** in 11.30 seconds. New tests compare local-window RGB SSIM against scikit-image, verify PSNR from MSE, check exact-image and dimension-mismatch behavior, and exercise the Streamlit image compression cards, three tabs and four chart outputs. Existing compression/framing and both FASTQ-source UI tests remain passing. The skipped test requires real ART. Production image metrics use SciPy, not a scikit-image runtime dependency.

# v13 single-page UI validation — 2026-10-06

The current revision exposes eight stages and HEDGES only. Streamlit AppTest checks both FASTQ sources, public settings import, consensus, decoding, and three-column summarization. Source switching does not reuse the other source's private file or reference.

Final full-suite run: **51 passed, 9 skipped** in 28.06 seconds. This includes the additional previous-step/summary test. UI-only suite: **4 passed, 1 skipped** in 8.20 seconds. Real ART is not installed and the DT4DDS import lacks Biopython in this execution environment; corresponding integration tests were skipped, not simulated using a replacement backend. The production DT4DDS and ART engines are unchanged.

The historical results below belong to earlier revisions, not new ART execution for v13.

# Historical DT4DDS + ART validation — 2026-09-29

## Result

**54 tests passed** in 35.63 seconds.

```bash
python -m pytest tests -q
```

Validation environment: Linux, Python 3.12, Streamlit 1.64.0, pytest 9.1.1, `dt4dds==1.1.0`, and real ART_Illumina Q Version 2.5.8 (June 6, 2016). ART was executed with its required GSL libraries; no fake FASTQ backend was used.

The Windows compatibility revision was rerun after replacing the evaluation-only `edlib` adapter with RapidFuzz and adding `install_windows.bat`: **54 tests passed in 38.73 seconds**. DT4DDS synthesis, aging and PCR code paths were unchanged.

The simplified Step 5.5 dashboard revision was rerun on 2026-09-29 UTC: **54 tests passed in 36.57 seconds**. Full-reference evaluation now uses RapidFuzz's compiled global edit operations; semi-global shortened-read evaluation retains the explicit dynamic-programming alignment. The decoder remains reference-blind.

The compact professor-facing simulation UI was rerun on 2026-09-30 UTC: **54 tests passed in 25.12 seconds**. AppTest exercised the equal-coverage selection through DT4DDS → real ART → consensus → exact HEDGES file recovery. Unit tests separately verify pool-proportional bias and equal per-strand allocation.

108 non-failing pandas migration warnings came from the retained original UI dataframe helper. They are unrelated to DT4DDS, ART, consensus or decoding.

## New workflow validation

- Official DT4DDS backend imports only version 1.1.0 and fails closed if unavailable.
- Twist synthesis, zero-exposure storage, Q5 PCR and finite physical sampling run deterministically for a fixed seed.
- Population abundance, unique variants, observed S/I/D burden and denominators are retained per stage.
- The manifest stores DT4DDS version, exact workflow settings and SHA-256 values for installed Twist, Aging, Q5 and properties profile files.
- DT4DDS molecular variants are passed to ART without resetting them to clean strands.
- ART executes in SE and PE modes, preserves sequence/quality lengths, produces matching mates, and is reproducible for a fixed seed.
- Variable-length molecular templates are grouped by length and retained; no input variant is silently discarded.
- FASTQ headers are neutral and do not contain SID or source-variant truth.
- ART SAM-derived sequencing S/I/D stays separate from molecular observed edit burden.
- Streamlit exposes one workflow, invalidates stale output after parameter changes and fails without leaving stale FASTQ when ART is missing.
- End-to-end AppTest recovered exact bytes `DNA test` through DT4DDS → real ART PE FASTQ → consensus → HEDGES decode.
- Page 2 package import, consensus and decoding remain operational.
- Existing compression, Hamming/no-Hamming SI, candidate search, search limits and reference-separation tests remain passing.

## Limits of this validation

This proves software execution and internal contracts, not calibration to a vendor or laboratory. Storage exposure remains zero in the locked experimental workflow. Only ART HS25 is integrated. Short-template handling is an explicit variable-read-length approximation, not adapter read-through. No Reed–Solomon or outer erasure code is present. Windows ART execution and browser screenshot layout were not tested in this Linux validation.
