# ART setup

ART is a separate executable, not a pip package. The app never fabricates FASTQ if ART is missing or fails.

## Windows / Anaconda

1. Download the Windows 64-bit ART package from the official distribution page:
   https://www.niehs.nih.gov/research/resources/software/biostatistics/art
2. Extract the entire package, keeping its companion files together.
3. In the Anaconda Prompt, test the binary (replace the example path):

```bat
"C:\tools\ART\art_illumina.exe" -h
```

4. Install the Windows-compatible DT4DDS core and start this app:

```bat
install_windows.bat
python -m streamlit run app.py
```

5. Step 5 → Advanced ART configuration → ART executable path: select/paste the actual executable path. Alternatively:

```bat
set ART_ILLUMINA=C:\tools\ART\art_illumina.exe
python -m streamlit run app.py
```

The Linux binary cannot run as a native Windows executable. If the native Windows distribution does not work on your machine, install ART and run **the whole Streamlit app in WSL/Linux**. This release does not translate Windows paths or dispatch individual commands through WSL/Docker. Windows native execution was not tested in this development environment.

## Linux / WSL

With the relevant distribution repositories available:

```bash
sudo apt install art-nextgen-simulation-tools
art_illumina -h
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Alternatively install an official Linux distribution and put `art_illumina` on PATH or set `ART_ILLUMINA` to its full path. Preserve any required runtime libraries. No ART binary is bundled with this project.

## What the app runs

Conceptual command for 20 pairs per 125-nt input amplicon:

```bash
art_illumina -ss HS25 -amp -p -sam -na -i input.fa -l 125 -c 20 -rs 12345 -qs 0 -qs2 0 -o reads_
```

`-c` is a count **per input sequence/amplicon**, not genomic fold coverage. For single-end, omit `-p` and `-qs2`. The wrapper groups templates by allocated count and effective read length and derives independent batch seeds from the displayed master seed. The exact commands and ART logs are available in Advanced and All results.

In the experimental workflow, sampling allocates the total read budget to DT4DDS sequence variants according to abundance BEFORE calling ART. Distinct variants are never flattened to equal coverage.

## ART sequencing indels on Windows

The default app setting is **Disabled — Windows-compatible default**. It passes `-ir 0 -dr 0` (and the R2 equivalents) to ART, while retaining ART's HS25 empirical quality profile and the associated substitution errors. It does not remove deletion/insertion mutations already present in the DT4DDS molecular template from synthesis or PCR.

Some Windows ART 2.5.8 builds abort with `qual size: ... read len: ...` when ART itself inserts or deletes a base in a read. Select **Use ART native defaults** only if your build has been independently validated. **Set custom rates** is provided for a validated build or controlled sensitivity analysis.

## Short templates and indels

ART can exit successfully while skipping a template shorter than the requested read length. This app validates both total and per-template counts and never accepts silently missing reads.

Two explicit choices are available:
- **Stop if template is short:** reject the run before invoking ART.
- **Shorten reads to template:** use `min(requested_length, template_length)` for that template. All sequence/quality characters still come from ART. This is a variable-read-length approximation, NOT a model of fixed-cycle adapter read-through or actual adapter trimming. The UI and metadata report how many read units were shortened. ART's indel handling near an amplicon boundary is retained, not overridden.

The original HEDGES oligos are 125 nt. Short SE reads may not span both primers and cannot necessarily be decoded by the current full-strand assignment pipeline. PE reads may be merged if overlap and primer constraints are met. Low recovery under an unsupported geometry must not be interpreted as HEDGES correction failure.

## Reproducibility and limits

Use the same ART binary, profile, master seed, inputs and app version. The run records executable SHA-256, exact batch commands, source-file hashes, template allocation, profile logs and configuration. Changing ART builds/platforms is not guaranteed to be bitwise identical.

Current interactive limits: 100,000 read units per run; nominal read length 25–150 nt, HS25 only; consensus uses at most 500 assigned units/SID. A capped-read warning appears if this limit is reached. Split larger studies into runs. Subprocess timeout, nonzero exit, malformed FASTQ, duplicate/unpaired reads and count mismatches fail explicitly. There is no fake-ART fallback.
