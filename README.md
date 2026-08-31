# ONT ITS Workflow - BRF @ Gadi

Pipeline for filtering, demultiplexing, orientation normalising and read-counting Oxford Nanopore ITS amplicon data on the NCI Gadi HPC.

---

## Overview

The workflow runs in a single PBS job (`run_script.qsub`) that calls two Python scripts in sequence:

| Stage | Script / Tool | What it does |
|-------|--------------|--------------|
| Pre-step | `generate_primer_setup.py` (run **twice**) | Converts the sample sheet CSV into two per-sample primer setup files -- one for Minibar, one for Cutadapt |
| Step 1 | `filter_chopper_demux_minibar.py` -> **Chopper** | Quality- and length-filters every raw `fastq.gz` |
| Step 2 | -> **Minibar** | Demultiplexes each filtered file into per-sample FASTQs, then merges across all `fastq_pass` files |
| Step 2.5 | -> **Cutadapt** (two-pass, `cutadapt_2pass.py`) | Normalises read orientation and trims 5'/3' constructs |
| Step 3 | -> organise + summarise | Groups reads by client, writes per-client and run-level read-count summaries |
| Final | -> log collection | Moves all `*.log` / `*.txt` files and `cutadapt_logs/` into `run_log_<date>/` |

After running the pipeline on the control sample and confirming that demultiplexing counts and per-sample yields match expectations, the same `run_script.qsub` can be submitted unchanged for production samples.

---

## Repository contents

```
run_script.qsub                  # PBS job wrapper -- edit Variables section here
generate_primer_setup.py         # Pre-step: sample sheet + barcode ref -> primer setup file
filter_chopper_demux_minibar.py  # Core pipeline: Steps 1, 2, 3, Final
cutadapt_2pass.py                # Step 2.5 two-pass cutadapt module (imported by core script)
Primer/                          # Barcode references + helpers (git-ignored, see below)
  IDT_ITS_96_barcode.txt         #   full construct -> cutadapt primer setup
  IDT_ITS_96_barcode_minibar.txt #   tail+index only -> Minibar primer setup
  make_minibar_reference.py      #   derives the *_minibar reference from the full one
  barcode_distance_qc.py         #   pairwise Hamming distance check on the index sets
```


---

## Prerequisites

| Tool | Default path |
|------|-------------|
| Python 3.12 | via `module load python3/3.12.1` |
| Chopper | `/g/data/vz35/ONT_16s_workflow/tools/chopper/chopper-linux-musl` |
| Minibar | `/g/data/vz35/ONT_16s_workflow/tools/minibar/minibar.py` |
| Cutadapt | `/g/data/vz35/ONT_16s_workflow/tools/cutadapt-env/bin/cutadapt` |
| Workflow scripts | `/g/data/vz35/ONT_ITS_workflow/tools/brf_ont_ITS_workflow/` |
| Barcode reference (cutadapt) | `/g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode.txt` |
| Barcode reference (Minibar) | `/g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode_minibar.txt` |

Chopper, Minibar and Cutadapt are shared with the 16S workflow and are **not** duplicated under `ONT_ITS_workflow/`; the paths above are the built-in defaults in the Python scripts and can be overridden with `--chopper`, `--minibar` and `--cutadapt`.

---

## Input files

### 1. Sample sheet (CSV)

| Column | Description |
|--------|-------------|
| `Client` | Client name -- used to organise output into subdirectories |
| `Sample_ID` | Sample identifier (alphanumeric, `-`, `_`; other characters are sanitised to `_`) |
| `Barcode` | Barcode ID matching a row in the barcode reference (`ITS01`-`ITS96`) |
| `Comment` | *(optional)* Shown in the per-client summary for samples below the low-read threshold |
| `Email` | *(optional)* Results sent through filesender -- not read by the pipeline |

Any sample whose `Barcode` is absent from the barcode reference aborts the pre-step with an error listing the offending rows.

### 2. Barcode references

Two TSVs describe the same 96 IDT ITS barcodes but differ in what the `FwPrimer` / `RvPrimer` columns contain:

| File | `FwPrimer` contains | Used by |
|------|--------------------|---------|
| `IDT_ITS_96_barcode.txt` | Illumina tail + 10 bp index + **ITS gene primer** (e.g. `...CGATCT` + `ATGACGTAGC` + `ACCWGCGGARGGATCATTA`) | Cutadapt -- the whole construct must come off the read |
| `IDT_ITS_96_barcode_minibar.txt` | Illumina tail + 10 bp index **only** (`...CGATCT` + `ATGACGTAGC`) | Minibar -- demultiplexing only needs to recognise the barcoded tail |

Both share the required columns `Barcode`, `FwIndex`, `FwPrimer`, `RvIndex`, `RvPrimer`; extra columns (`Tem_index_F`, `Tem_index_R`) are ignored.

Splitting the reference this way keeps the IUPAC-degenerate gene primer (`W`, `R`, `S`, `N`, `D`) out of Minibar's matching, since Minibar scores each degenerate position against a single concrete base and would otherwise spend the mismatch budget on primer design rather than read quality.

During later test, this seems not influence the minibar function but this change still keeps.

The `*_minibar` file is **derived, not hand-maintained**: `Primer/make_minibar_reference.py`
cuts `FwPrimer` / `RvPrimer` at the end of the index and carries every other column through
unchanged.

```bash
python3 Primer/make_minibar_reference.py Primer/IDT_ITS_96_barcode.txt
# -> Primer/IDT_ITS_96_barcode_minibar.txt
```

Re-run it whenever the full-construct reference changes, so the two stay in step.

### 3. Raw reads

PromethION `fastq_pass` directory located under:

```
/g/data/vz35/PromethION_data/sequencer_uploads/<run_name>/
```

The wrapper locates `fastq_pass` automatically from `run_name` and aborts if zero, or more than one, are found.

---

## Configuration

Edit the **Variables** section near the top of `run_script.qsub`:

```bash
# --- Variables (edit per run) ---
run_name=ONT_ITS_20260617          # Subdirectory under sequencer_uploads/
samplesheet=/g/data/vz35/ONT_ITS_workflow/sample_sheet/ITS_samplesheet_20260721.csv

# Same barcode set in two forms, one per tool
barcode_file_minibar=/g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode_minibar.txt
barcode_file_cutadapt=/g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode.txt

output_dir=/g/data/vz35/ONT_ITS_workflow/minibar_output/ONT_ITS_TBC_${date_stamp}

# --- Variables index, change as needed ---
min_quality=15                     # Chopper -q
min_length=400                     # Chopper --minlength (ITS default)
max_length=1200                    # Chopper --maxlength (ITS default)

primer_mismatch=9                  # Minibar -E
barcode_mismatch=1                 # Minibar -e
```

**Both** barcode references are set explicitly and passed to `generate_primer_setup.py` with `-b`, so nothing depends on the generator's built-in default. Both files must cover every barcode ID used in the sample sheet's `Barcode` column.

If you swap in a different barcode set, also re-check `min_length` / `max_length` and `primer_mismatch` -- both are tuned for the IDT ITS construct.

Everything in the `DONT-CHANGE` section below resolves automatically from the variables above: `raw_input_dir` is found from `run_name`, and the two generated primer setup files are derived from `primer_output_dir` + `date_stamp`.

---

## Usage

### 1. Prepare the sample sheet

Fill in the CSV with `Client`, `Sample_ID` and `Barcode` columns.
Include a control sample with a known expected yield to validate the run.

Optional `Comment` can be put based on the amplification tendency.

### 2. Update the Variables section

Open `run_script.qsub` and set `run_name`, `samplesheet` and `output_dir`.
If the barcode set has changed, point `barcode_file_minibar` and `barcode_file_cutadapt` at the new pair and re-check `min_length` / `max_length` and `primer_mismatch`.

### 3. Submit the job

```bash
qsub run_script.qsub
```

PBS resources requested: 2 CPUs, 10 GB RAM, 10 GB jobfs, 20 h walltime, project `vz35`, queue `biodev`, storage `gdata/vz35`.

### 4. Validate with the control run

After the job completes, open `run_log_<date>/read_counts_summary.txt` and the control sample's client `summary.txt` inside `integrated_demultiplexing/<Client>/`.
Confirm that the control sample's read count and percentage match the expected values for that barcode, and check `run_log_<date>/cutadapt_summary.txt` for a plausible kept / discarded / reverse-complemented split.
Once validated, the same script can be resubmitted for any other sample set by updating the Variables section.

### 5. (Optional) Generate the primer setup files standalone

Both files are generated automatically inside the job. To create them independently:

```bash
# Minibar file (tail + index only) -- built first, then renamed out of the way
python3 generate_primer_setup.py <samplesheet.csv> \
    -b /g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode_minibar.txt \
    -o OUTPUT_DIR --date YYYYMMDD
mv OUTPUT_DIR/ITS_primer_setup_YYYYMMDD.txt OUTPUT_DIR/ITS_primer_setup_minibar_YYYYMMDD.txt

# Cutadapt file (full construct) -- keeps the default name
python3 generate_primer_setup.py <samplesheet.csv> \
    -b /g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode.txt \
    -o OUTPUT_DIR --date YYYYMMDD
```

This is exactly what the wrapper does: the same generator invoked twice against different `-b` inputs. The generator always writes `ITS_primer_setup_<date>.txt`, which is why the Minibar file is built first and renamed before the cutadapt file is built. Both are tab-separated with columns:
`SampleID`, `FwIndex`, `FwPrimer`, `RvIndex`, `RvPrimer`.

`-b` defaults to `IDT_ITS_96_barcode.txt`, but the wrapper always passes it explicitly -- supply a different file for runs with external barcodes.

---

## Pipeline steps in detail

### Pre-step: generate_primer_setup.py

Reads the sample sheet and looks up each barcode in the supplied reference to produce a tab-separated primer setup file. Run twice by the wrapper to produce:

| File | Built from (`-b`) | Consumed by |
|------|------------------|-------------|
| `ITS_primer_setup_minibar_<date>.txt` | `barcode_file_minibar` | Step 2 (Minibar) |
| `ITS_primer_setup_<date>.txt` | `barcode_file_cutadapt` | Step 2.5 (Cutadapt) |

Both are written to `sample_sheet/` and passed to the core script as `--minibar-primer-file` and `--primer-file` respectively. If `--minibar-primer-file` is omitted, Minibar falls back to `--primer-file`.

---

### Step 1: Chopper -- quality and length filtering

Each `*.fastq.gz` file in `fastq_pass/` is piped through Chopper, and reads that pass all three filters are retained:

| Filter | Parameter | Value |
|--------|-----------|-------|
| Quality score | `-q` | >= 15 |
| Minimum length | `--minlength` | 400 bp |
| Maximum length | `--maxlength` | 1,200 bp |

All three are set in the **Variables** section of `run_script.qsub` (`min_quality`, `min_length`, `max_length`) and passed through to the core script.

The length window must bracket the expected amplicon for the primer pair in use, so re-check it whenever the primers change. The 400-1,200 bp default is sized for ITS: the main peak sits at ~0.78-0.80 kb with a shoulder at ~0.60-0.70 kb, and concatemers/junk appear above ~1.5 kb.

The equivalent command run per file:

```
zcat <file>.fastq.gz | chopper -q 15 --minlength 400 --maxlength 1200 | gzip > <file>_filtered.fastq.gz
```

Filtered files are written to `chopper_filtered/` and deleted automatically after demultiplexing.

---

### Step 2: Minibar -- demultiplexing

Minibar demultiplexes each filtered file into per-sample FASTQs using the **Minibar** primer setup file.
Each file is processed in its own subdirectory, and per-sample results are merged across all `fastq_pass` files at the end.

| Flag | Value | Set in | Meaning |
|------|-------|--------|---------|
| `-e` | 1 | `barcode_mismatch` (qsub) | Allowed edit distance on the barcode index |
| `-E` | 9 | `primer_mismatch` (qsub) | Allowed edit distance on the primer |
| `-l` | 200 | hard-coded in Python | Search window at each read end (bp) |
| `-M` | 2 | hard-coded in Python | Require barcode match on both ends |
| `-F` | -- | hard-coded in Python | Write each sample to its own file |

**`-M 2` must not be relaxed.** The IDT ITS 96 barcode set is *combinatorial*, not unique-dual: 12 forward indices x 8 reverse indices generate the 96 samples, so each `FwIndex` is shared by 8 samples and each `RvIndex` by 12. Only the *pair* identifies a sample. Under `-M 1`, a single-end match would be ambiguous across 8-12 candidates and reads would be silently cross-assigned between samples instead of landing in `sample_unk`.

**`-e` (barcode mismatch) is kept at 1.** Raising it above 1 can rescue some unassigned (`unk`) reads but also grows `Multiple_Matches`, because widening tolerance can only add candidate matches, never remove them, and the IDT ITS `FwIndex` set has a minimum pairwise Hamming distance of only 3 (`RvIndex` is 5). If it is ever raised, watch whether the `Multiple_Matches` growth concentrates on the known close pairs: F1/F11, F1/F10, F3/F9, F6/F8, F8/F10, F9/F12, F10/F11 (all distance <= 4). Use `Primer/barcode_distance_qc.py` to re-derive these distances for a new barcode set.
*Tested 20260722 -- keep `barcode_mismatch=1`.*

**`-E` (primer mismatch) is set to 9 for ITS.** The value is a hangover from matching against the full degenerate construct, where every IUPAC position counted as a mismatch before any real sequencing error. Now that Minibar reads the tail-and-index-only reference, the constant Illumina tail carries no degeneracy and the budget is less critical.
*Tested 20260722 -- `primer_mismatch` 5 and 9 gave the same result.* Re-check whenever the primers or the Minibar barcode reference change.

`-l 200` is the per-end search window. The barcoded construct is ~32 bp in the Minibar reference (22 bp Illumina tail + 10 bp index), so 200 leaves ample headroom for the few bases of noise ONT reads often carry before the adapter begins.

Output files: `sample_<SampleID>.fastq`, `sample_unk.fastq`, `sample_Multiple_Matches.fastq`.
Per-file Minibar subdirectories are removed after merging to save storage.

---

### Step 2.5: Cutadapt -- orientation normalisation and primer trimming

Run on each `sample_<SampleID>.fastq` after merging, using the **cutadapt** primer setup file (full construct). Reads from `sample_unk` and `sample_Multiple_Matches` are left untouched, as are samples with no matching row in the primer setup file.

**Pass 1 -- orient and trim 5' construct (strict)**

```
cutadapt -g FWD_PRIMER --revcomp --rename={header} --discard-untrimmed -e 0.2 -j 4 --json ...
```

- Searches for the forward construct at the 5' end `-g FWD`
- If it is found on the minus strand, the read is reverse-complemented `--revcomp` so all output reads face the same direction
- `--rename={header}` keeps the read name unchanged (cutadapt otherwise appends ` rc` to reverse-complemented reads)
- Reads where the forward construct cannot be found are discarded `--discard-untrimmed` -- these are likely noise or off-target sequences

**Pass 2 -- trim 3' construct (tolerant)**

```
cutadapt -a REV_PRIMER_RC -e 0.2 -j 4 --json ...
```

- Trims the reverse-complemented reverse construct from the 3' end (reverse complement is IUPAC-aware)
- Reads where the 3' construct is not found are **kept** (truncated reads are retained)

Output: `<SampleID>.fastq.gz` replaces `sample_<SampleID>.fastq`.
Per-sample logs go to `cutadapt_logs/`; the JSON reports are parsed for stats and then deleted. A run-level `cutadapt_summary.txt` reports per-sample input / kept / discarded / reverse-complemented / 3'-trimmed counts and percentages. A sample that fails cutadapt is reported as `ERROR` in the summary and its untrimmed `sample_*.fastq` is left in place.

---

### Step 3: Organise by client and summarise

Demultiplexed files are moved into `integrated_demultiplexing/<Client>/` based on the sample sheet. If cutadapt skipped a sample, the pre-cutadapt `sample_<id>.fastq` is moved instead and flagged in the log as `Moved (uncut)`.

Two summary files are written:

| File | Contents |
|------|----------|
| `integrated_demultiplexing/<Client>/summary.txt` | Per-sample read counts for that client; shows the `Comment` column for samples below 15,000 reads |
| `read_counts_summary.txt` | Run-level totals: total filtered input, total demultiplexed, successfully demultiplexed, and per-sample percentages |

`SUCCESSFULLY DEMULTIPLEXED` is computed as total filtered input minus `sample_unk` minus `sample_Multiple_Matches`.

The low-read threshold is the `LOW_READ_THRESHOLD` constant in `filter_chopper_demux_minibar.py` (set to `None` to disable).

### Final: log collection

All top-level `*.log` and `*.txt` files (including `read_counts_summary.txt`, `cutadapt_summary.txt` and the `cutadapt_logs/` directory) are moved into `run_log_<date>/` to keep the output root clean. Per-client `summary.txt` files stay with their client directory.

---

## Output structure

```
minibar_output/ONT_ITS_TBC_<date>/
├── integrated_demultiplexing/
│   ├── <ClientA>/
│   │   ├── <SampleID>.fastq.gz
│   │   └── summary.txt
│   ├── <ClientB>/
│   │   └── ...
│   ├── sample_Multiple_Matches.fastq
│   └── sample_unk.fastq
└── run_log_<date>/
    ├── read_counts_summary.txt
    ├── cutadapt_summary.txt
    └── cutadapt_logs/
        ├── <SampleID>.pass1_orient.log
        └── <SampleID>.pass2_trim3.log
```

- `integrated_demultiplexing/` -- reads merged across all `fastq_pass` files and split by client
- `sample_unk.fastq` -- reads that did not match any barcode
- `sample_Multiple_Matches.fastq` -- reads that matched more than one barcode
- `summary.txt` -- per-client read counts; low-count samples (< 15,000 reads) show the `Comment` column from the sample sheet
- `read_counts_summary.txt` -- overall run summary with percentages relative to total filtered input
- `chopper_filtered/` -- intermediate, deleted after demultiplexing

---

## Tool parameter reference

### Chopper

| Parameter | Value | Effect |
|-----------|-------|--------|
| `-q 15` | 15 | Minimum mean quality score |
| `--minlength 400` | 400 bp | Discard reads shorter than this |
| `--maxlength 1200` | 1,200 bp | Discard reads longer than this |

### Minibar

| Parameter | Value | Effect |
|-----------|-------|--------|
| `-e 1` | 1 | Barcode index edit distance -- keep at 1, see Step 2 |
| `-E 9` | 9 | Primer edit distance (script default 5; wrapper passes 9 for ITS) |
| `-l 200` | 200 | End-window search length (bp) |
| `-M 2` | 2 | Match barcodes on both ends -- required, see Step 2 |
| `-F` | -- | Write separate file per sample |

### Cutadapt

| Parameter | Value | Effect |
|-----------|-------|--------|
| `-e 0.2` | 0.2 | Error rate for adapter matching (both passes) |
| `-j 4` | 4 | Threads per sample |
| `--revcomp` | -- | Pass 1: reverse-complement reads where the construct is found on the minus strand |
| `--rename={header}` | -- | Pass 1: preserve the original read name |
| `--discard-untrimmed` | -- | Pass 1: drop reads with no forward construct match |
| `-g FWD` | forward construct | Pass 1: 5' adapter |
| `-a REV_RC` | RC of reverse construct | Pass 2: 3' adapter |

---

## Core script CLI

`run_script.qsub` calls the core script with the options below; they can also be used directly for a manual re-run.

| Option | Default | Notes |
|--------|---------|-------|
| `--raw-input-dir` | *required* | Directory of raw `*.fastq.gz` (i.e. `fastq_pass`) |
| `--output-dir` | *required* | Run output directory (created if absent) |
| `--primer-file` | *required* | Cutadapt primer setup file |
| `--minibar-primer-file` | falls back to `--primer-file` | Minibar primer setup file |
| `--samplesheet` | *required* | Sample sheet CSV |
| `--chopper` | Gadi 16S tools path | Chopper binary |
| `--minibar` | Gadi 16S tools path | `minibar.py` |
| `--cutadapt` | Gadi 16S tools path | Cutadapt binary |
| `--min-quality` | 15 | Chopper `-q` |
| `--min-length` | 400 | Chopper `--minlength` |
| `--max-length` | 1200 | Chopper `--maxlength` |
| `--barcode-mismatch` | 1 | Minibar `-e` |
| `--primer-mismatch` | 5 | Minibar `-E` (wrapper passes 9) |
| `--cutadapt-threads` | 4 | Cutadapt `-j` |
| `--cutadapt-error-rate` | 0.2 | Cutadapt `-e` |
