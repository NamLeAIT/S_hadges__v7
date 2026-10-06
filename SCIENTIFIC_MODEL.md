# Scientific model and claim boundary

## Primary sources

- Gimpel AL, Stark WJ, Heckel R, Grass RN. *A digital twin for DNA data storage based on comprehensive quantification of errors and biases*. Nature Communications 14, 6026 (2023). DOI: https://doi.org/10.1038/s41467-023-41729-1
- DT4DDS author implementation: https://github.com/fml-ethz/dt4dds
- Huang W et al. *ART: a next-generation sequencing read simulator*. Bioinformatics 28, 593–594 (2012). DOI: https://doi.org/10.1093/bioinformatics/btr708

The runtime pins `dt4dds==1.1.0`. Each run stores the package version and SHA-256 values of the four DT4DDS profile files used by the workflow.

## Mathematical workflow

For designed central oligo `s_i`, DT4DDS synthesis draws a sequence abundance using the official Twist lognormal profile. A finite pool sample is then mutated by the official synthesis error generators. The app does not replace the source nucleotide/repeat biases with independent uniform errors.

Primers are attached according to the DT4DDS reference workflow. Storage uses `Aging`. At exposure `h=0`, the pool passes unchanged; for a future nonzero exposure, the source model scales its error rates by half-lives and applies its base-specific spectrum. Half-life is a normalized exposure coordinate, not hydrogen, years or temperature.

Q5 PCR uses DT4DDS's cycle-level binomial amplification counts and sequence-specific efficiency property. Its aggregate substitution generator uses:

`mean rate = cycles × polymerase base-substitution rate / polymerase fidelity`

with the official Q5 substitution spectrum. A mutation already present in a parent sequence remains in amplified descendants. The simulator may compress high copy-number populations internally while preserving represented counts; therefore a represented molecule is not always one explicitly materialized Python object.

After PCR, the app draws a finite abundance-weighted aliquot. ART read units are allocated from that aliquot in proportion to represented molecule count. ART receives each actual mutated sequence, including variable lengths caused by molecular indels.

ART is run with `-ss HS25 -amp -sam -na`. Input variants are grouped only by requested count and effective read length because ART accepts one read length per invocation. The wrapper verifies output counts, mate names, sequence/quality length equality and SAM mappings. It then replaces headers with neutral IDs and shuffles read units; it does not rewrite bases or Phred scores.

## Stage measurement

DT4DDS returns population states, not the causal mutation history of every represented descendant. Therefore stage molecular error counts are computed as follows:

1. assign each observed unique sequence to its nearest designed strand using global edit distance;
2. retain deterministic tie handling and report abundance assigned under ties;
3. compute S/I/D from a minimum-edit path;
4. multiply each edit count by represented molecular abundance;
5. divide by designed reference nucleotides × represented molecule count.

These are **observed edit-burden measurements**, not new-event counts. They can decrease after sampling and can multiply after PCR because a mutated lineage has many descendants. Stage denominators are recorded.

ART sequencing errors are different: they are measured from ART's SAM CIGAR against the exact input molecular template and are new sequencing errors. Molecular and sequencing error counts must not be blindly added because their denominators and meanings differ.

## Decoder separation

- SID assignment used for plots and abundance tables is evaluation-only.
- FASTQ names do not contain SID, source variant or clean reference.
- Consensus receives only R1/R2 sequence, Phred quality and explicit public design settings.
- HEDGES candidate search receives only consensus/output settings.
- Reference comparisons occur after processing.
- Exact recovery uses SHA-256 of the stored representation and file reconstruction.

This separation prevents simulated ground truth from becoming a decoder shortcut.

## Fixed workflow versus calibration

The UI locks the molecular profile to make demonstrations understandable and comparable. Locked does not mean universally correct. Current workflow choices include 100 synthesis-sampled molecules/strand, 0 storage half-lives, Q5 PCR with 20 cycles and 0.95 mean efficiency, and a 100-molecule/strand physical aliquot budget.

For a lab-specific model, reserve independent experimental data and estimate at least:

- synthesis abundance distribution and context-dependent indels;
- pre-/post-storage damage under controlled time, temperature and humidity;
- PCR abundance variance across low/medium/high cycles with technical replicates;
- positional R1/R2 Phred distributions and substitution/indel spectra;
- zero-coverage and low-coverage SID probability after physical sampling.

Fit on a calibration subset and validate on held-out FASTQ. Report parameter uncertainty and sensitivity, not only a single best fit.

## Claims this build supports

- a real execution of DT4DDS 1.1.0 molecular process models followed by real ART execution;
- a reproducible configuration and provenance manifest;
- inherited/cascaded molecular errors and abundance-weighted sampling;
- FASTQ with real ART bases and Phred qualities;
- reference-free consensus/HEDGES processing followed by evaluation;
- exact-file recovery experiments over seeds and coverage.

## Claims this build does not support

- calibration to a particular vendor, Korean synthesis company, PCR kit or Illumina instrument;
- exact mutation genealogy per molecule;
- conversion from half-lives to calendar time without an independently measured decay law;
- correction of a completely missing strand by inner HEDGES alone;
- a formal arbitrary-K HEDGES guarantee;
- 99% recovery from a small number of Monte Carlo trials;
- research novelty solely from integrating DT4DDS, ART or Streamlit.
