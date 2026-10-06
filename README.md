# DNA Data Storage Simulator — DT4DDS + ART

Two-page Streamlit application for an end-to-end DNA data-storage experiment:

`file → compression → DNA/HEDGES encoding → DT4DDS molecular channel → ART Illumina FASTQ → consensus → HEDGES decode → file verification`

No Reed–Solomon or outer erasure code is included. The app does not silently replace DT4DDS or ART with a homemade fallback.

## Install and run

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

On native Windows, use the supplied installer from an activated Python 3.10 environment:

```bat
install_windows.bat
python -m streamlit run app.py
```

The Windows installer deliberately skips DT4DDS's optional `pysam` and `edlib` analysis dependencies because their current PyPI distributions do not provide the required native Windows wheels. The app uses DT4DDS's synthesis, aging and PCR processes, none of which import those analysis modules. Evaluation alignment uses RapidFuzz's native Windows wheel. This changes installation only; it does not replace or reimplement the DT4DDS molecular processes.

ART is a separate executable. Install it using `ART_SETUP.md`, put `art_illumina` on `PATH`, select the executable in Advanced settings, or set `ART_ILLUMINA`.

## Compression and file framing

The compression stage benchmarks the uploaded file together with domain-aware candidates:

- text/structured data: Brotli, Zstandard and XZ;
- general binary, PDF and Office data: Zstandard, XZ and zlib;
- image: optimized PNG, lossless WebP, and JPEG/WebP at Light, Medium and Strong levels;
- audio: FLAC and Opus at 128/64/32 kbit/s;
- video: H.264 and H.265 at Light, Medium and Strong CRF levels.

Image candidates report PSNR/SSIM/MAE, audio candidates report decoded-waveform SNR/correlation/NRMSE, and video candidates report FFmpeg PSNR/SSIM. The UI previews the original and selected image, audio, video, text, PDF or extractable OOXML content. General lossless candidates use a self-describing `DDSC` container carrying the original name, MIME type, byte length and CRC32.

### Compression quality dashboard (v13.1)

The selected output shows PSNR, SSIM and MAE cards for images, or the corresponding audio/video metrics. Three comparison tabs show quality per compressor, stored size versus measured quality, and file sizes. Green bars mark the selected candidate; exact (infinite dB) results remain in the table and cards but are omitted from the dB plot. SSIM × 100 is a similarity score, not file-recovery probability.

Image SSIM uses local uniform windows up to 7 × 7, RGB channel averaging, sample covariance and data range 255, consistent with the default scikit-image local-window definition. Dimension mismatches are not resized silently. Tiny images below 3 pixels on either axis have no SSIM. Text/PDF/binary display file sizes rather than inapplicable image metrics. Compression quality is separate from later DNA decoding integrity.

The interface exposes HEDGES only. Simple Mapping and rule-based encoders are not selectable. JPEG and WebP use quality 90/70/10 (Light/Medium/Strong), matching the supplied compressor file; the additional self-describing exact compressors are retained.

New HEDGES files embed the stored file name, domain, exact byte length and CRC32 inside the encoded bitstream before 103-bit fragmentation. Padding removal and file reconstruction therefore do not depend on Streamlit session state or a hidden copy of the original file.

## Single-page workflow

The eight sections are Input, Compression, Encoding, NGS Strand Preparation, FASTQ Generation, Consensus, Decoding and Summarization. NGS means next-generation sequencing. Step 5 retains the experimental workflow:

1. **Twist array synthesis — DT4DDS 1.1.0**
   - Official `ArraySynthesis_Twist` profile.
   - Lognormal abundance sigma 0.30.
   - 100 sampled molecules per designed strand before primer attachment.
   - The profile includes its measured deletion rate, nucleotide bias and repeat-context bias.
2. **Storage — DT4DDS Aging**
   - Default exposure is 0 half-lives because the current laboratory plan is synthesis → PCR → sequencing, not an aging experiment.
   - The stage remains in the pipeline and manifest; it is not silently removed.
3. **PCR — DT4DDS Q5**
   - 20 cycles, mean efficiency 0.95, template/reaction volume 1/20 µL.
   - Sequence-specific amplification efficiency and inherited molecular errors are retained by DT4DDS.
4. **Physical sampling**
   - A 100-molecule-per-designed-strand aliquot is sampled from the PCR pool by abundance.
   - Low-abundance SIDs may be absent; no artificial minimum coverage is imposed.
5. **Illumina sequencing — ART HS25**
   - The user selects mean read coverage, single/paired-end and read length.
   - ART generates FASTQ bases and Phred qualities. The app preserves them and verifies them against ART SAM.
   - ART-generated sequencing indels are disabled by default for compatibility with affected Windows ART 2.5.8 builds; ART's empirical Q-derived substitutions remain enabled. Molecular indels already present in DT4DDS templates are retained.

The fixed molecular workflow is intentionally simple in the UI. Exact parameters, package version, configuration-file SHA-256 values, ART executable SHA-256 and ART commands are saved in `channel_manifest.json`.

## What the measurements mean

- **Simulation result:** sampled-strand count, missing-strand count, actual mean coverage and total FASTQ edit rate before consensus.
- **Per-strand table:** coverage, substitution, insertion, deletion, total edits and error rate for every designed strand. This is the only main result table; process-stage error charts are retained only in downloaded audit files.
- **Sampling choice:** pool-distribution sampling follows molecular abundance; equal-coverage sampling assigns the same target read units to each surviving strand while retaining abundance-weighted molecular variants within that strand.
- **Molecular stage S/I/D:** abundance-weighted minimum-edit burden against the nearest designed strand. DT4DDS does not expose a per-molecule mutation genealogy, so these values are not labeled “new events.” Inherited PCR errors may be counted in many descendant molecules.
- **ART sequencing S/I/D:** new sequencing errors measured directly from ART's explicit `=/X/I/D` SAM CIGAR against the sampled molecular template.
- **Total read error:** post-hoc alignment of observed reads to the clean design for evaluation only.
- **Consensus improvement:** BODY edit rate before/after consensus over the same nonempty SID groups; missing groups are reported separately.
- **Final success:** exact stored-byte SHA-256 and reconstructed-file checks.

The clean reference, DT4DDS-to-SID evaluation map and ART SAM do not enter consensus or HEDGES candidate selection. FASTQ headers are neutral (`read000...`) and contain no SID.

## Consensus input

Step 6 places the previous-step FASTQ source on the left and external R1/R2 FASTQ upload on the right. Select one source, explicitly apply any uploaded public settings, then run:

`FASTQ → assignment/alignment → consensus → HEDGES search → file reconstruction`

HEDGES keeps `OK`, `FAIL`, `AMBIGUOUS` and `SEARCH_LIMIT` separate. Increasing candidate/time budgets may resolve `SEARCH_LIMIT`, but cannot turn genuine ambiguity into a guarantee.

Step 8 shows Original / Compressed / Decoded previews, selected compression metrics, strand properties and recovery integrity. External FASTQ does not inherit the current input file or its private reference. Missing reference measurements remain unavailable. Lossy compression is evaluated separately from exact stored-byte recovery.

## Downloads

The single **Download FASTQ package** button includes:

- `R1.fastq` and optional `R2.fastq`;
- `decoder_settings.json`;
- `channel_manifest.json`;
- evaluation-only reference, expected hash and configuration files.

Secondary audit artifacts remain under **All results**.

## Reproducible experiments

```bash
python run_channel_experiment.py \
  --strands prepared_strands.csv \
  --settings decoder_settings.json \
  --expected-hash STORED_SHA256 \
  --mode experimental \
  --coverage 5,10,20,50,100,350 \
  --repeats 30 \
  --paired \
  --output experiment_output
```

The runner reports exact-file recovery, operational failures, consensus errors and Wilson 95% confidence intervals. A few successful trials do not establish 99% reliability.

## Scientific boundary

This build is a reproducible **published-profile simulation**, not yet a calibrated digital twin of a particular synthesis company, PCR protocol, Illumina run or laboratory. Use held-out synthetic-DNA FASTQ and process metadata to fit and validate abundance, substitution/indel and positional-quality distributions before making lab-specific claims. See `SCIENTIFIC_MODEL.md`.

## Tests

```bash
python -m pytest tests -q
```

Real ART tests are skipped when ART is unavailable. `tests/test_dt4dds_backend.py` verifies the official DT4DDS backend, reproducibility, abundance and manifest provenance.

## Dependency license note

DT4DDS 1.1.0 is distributed under GPLv3. This project references it as an external pinned dependency rather than copying its source. Review the upstream license and your distribution obligations before publishing or redistributing a packaged deployment.
