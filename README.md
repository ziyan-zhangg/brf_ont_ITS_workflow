# ONT ITS Workflow - BRF @ Gadi

Pipeline for filtering, demultiplexing, orientation normalising, read-counting and per-sample QC of Oxford Nanopore ITS amplicon data on the NCI Gadi HPC.

---

## Overview

The workflow runs in a single PBS job (`run_script.qsub`) that calls a pre-step script, the core pipeline and the QC report:

| Stage | Script / Tool | What it does |
|-------|--------------|--------------|
| Pre-step | `generate_primer_setup.py` | Converts the sample sheet CSV into two per-sample primer setup files -- one for Minibar, one for Cutadapt |
| Step 1 | `filter_chopper_demux_minibar.py` -> **Chopper** | Quality- and length-filters every raw `fastq.gz` — **`pre` mode only** |
| Step 2 | -> **Minibar** | Demultiplexes each input file into per-sample FASTQs, then merges across all `fastq_pass` files |
| Step 3 | -> `demux_qc.py` | Per-sample length/quality profile of the demultiplexed reads, before filtering |
| Step 4 | -> **Chopper** | Quality- and length-filters each merged per-sample FASTQ — **`post` mode only** |
| Step 5 | -> **Cutadapt** (two-pass, `cutadapt_2pass.py`) | Normalises read orientation and trims 5'/3' constructs |
| Step 6 | -> organise + summarise | Groups reads by client, writes per-client and run-level read-count summaries |
| Step 7 | -> log collection | Moves all `*.log` / `*.txt` files, `cutadapt_logs/` and `minibar_logs/` into `run_log_<date>/` |
| Post-step | `qc_report.py` | Per-sample PASS / CHECK / FAIL report into `qc_report/` (see [Per-sample QC report](#per-sample-qc-report)) |

After running the pipeline on the control sample and confirming that demultiplexing counts and per-sample yields match expectations, the same `run_script.qsub` can be submitted unchanged for production samples.

### Filter stage: when Chopper runs

`filter_stage` in `run_script.qsub` decides whether Chopper runs **before** demultiplexing (step 1) or **after** it, per sample (step 4):

| Mode | Chopper runs | Minibar sees | Notes |
|------|--------------|--------------|-------|
| `post` *(default)* | Step 4, on merged per-sample FASTQs | Raw `fastq_pass` reads | No pre-filter at all — `fastq_pass` is already MinKNOW filtered |
| `pre` | Step 1, on whole `fastq.gz` files | Filtered reads | Legacy ordering (the only one before this release) |

Both modes share one step numbering, so step *N* means the same thing in either job log and the two can be read side by side. Each mode skips exactly one step, and the skip is announced in the log rather than silently absent:

```
pre                                   post
STEP 1: Chopper filtering             STEP 1: Chopper filtering SKIPPED
STEP 2: Minibar demultiplexing        STEP 2: Minibar demultiplexing
STEP 3: Per-sample demux QC           STEP 3: Per-sample demux QC
STEP 4: Per-sample Chopper SKIPPED    STEP 4: Chopper filtering (per sample)
STEP 5: Cutadapt two-pass             STEP 5: Cutadapt two-pass
STEP 6: Organising reads by client    STEP 6: Organising reads by client
STEP 7: Collecting logs               STEP 7: Collecting logs
```

**Why `post` is the default.** Filtering before demultiplexing throws reads away before anyone can see what they looked like, so a low per-sample yield is ambiguous: too few reads, or the right reads outside the assumed 400–1,200 bp window? Running Minibar first makes the demux loss and the filter loss separately visible (steps 3 and 4).

On ITS runs MinKNOW applies a **400 bp** minimum read length and pass Q15, the same floor as Chopper (`-q 15 --minlength 400 --maxlength 1200`), so in practice Chopper only removes reads > 1,200 bp. Post mode therefore adds only those > 1.2 kb reads to Minibar's load, and in exchange the QC report can measure them per sample (`long_frac`). `-M 2` still applies to them, so they are only assigned on a both-ends barcode match.

Both `run_script.qsub` and the core script default to `post`.

---

## Repository contents

```
run_script.qsub                  # PBS job wrapper -- edit Variables section here
generate_primer_setup.py         # Pre-step: sample sheet + barcode ref -> primer setup file
filter_chopper_demux_minibar.py  # Core pipeline: Steps 1-7
demux_qc.py                      # Step 3: per-sample demux QC (imported by core script)
cutadapt_2pass.py                # Step 5: two-pass cutadapt (imported by core script)
qc_report.py                     # Post-step: per-sample QC report (reads demux_qc/ + counts)
qc_report_template.html          # Page template for qc_report.html (must sit beside qc_report.py)
qc_config.yaml                   # QC report thresholds for ITS, versioned
run_qc_report.qsub               # PBS wrapper to re-run qc_report.py on a finished run
tests/test_qc_report.py          # Fixture test for qc_report.py / demux_qc.py
Primer/                          # Barcode references + helpers (git-ignored, see below)
  IDT_ITS_96_barcode.txt         #   full construct -> cutadapt primer setup
  IDT_ITS_96_barcode_minibar.txt #   tail+index only -> Minibar primer setup
  make_minibar_reference.py      #   derives the *_minibar reference from the full one
  barcode_distance_qc.py         #   pairwise Hamming distance check on the index sets
```

`demux_qc.py` and `cutadapt_2pass.py` are imported as siblings of the core script, and `qc_report.py` imports the core script, so all five `.py` files plus `qc_report_template.html` and `qc_config.yaml` must sit in the same directory on Gadi.

`qc_report.py`, `qc_report_template.html` and `demux_qc.py` are shared with the 16S workflow. Everything assay-specific lives in `qc_config.yaml`, so the 16S repo can take the same `qc_report.py` / template unchanged.

---

## Prerequisites

| Tool | Default path |
|------|-------------|
| Python 3.12 | via `module load python3/3.12.1` |
| Chopper | `/g/data/vz35/ONT_16s_workflow/tools/chopper/chopper-linux-musl` |
| Minibar | `/g/data/vz35/ONT_16s_workflow/tools/minibar/minibar.py` |
| Cutadapt | `/g/data/vz35/ONT_16s_workflow/tools/cutadapt-env/bin/cutadapt` |
| QC venv (PyYAML, matplotlib) | `/g/data/vz35/ONT_16s_workflow/tools/qc_env` |
| Workflow scripts | `/g/data/vz35/ONT_ITS_workflow/tools/brf_ont_ITS_workflow/` |
| Barcode reference (cutadapt) | `/g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode.txt` |
| Barcode reference (Minibar) | `/g/data/vz35/ONT_ITS_workflow/tools/IDT_ITS_96_barcode_minibar.txt` |

Chopper, Minibar, Cutadapt and the QC venv are shared with the 16S workflow and are **not** duplicated under `ONT_ITS_workflow/`; the tool paths above are the built-in defaults in the Python scripts and can be overridden with `--chopper`, `--minibar` and `--cutadapt`.

---

## Input files

### 1. Sample sheet (CSV)

| Column | Description |
|--------|-------------|
| `Client` | Client name -- used to organise output into subdirectories |
| `Sample_ID` | Sample identifier (alphanumeric, `-`, `_`; other characters are sanitised to `_`) |
| `Barcode` | Barcode ID matching a row in the barcode reference (`ITS01`-`ITS96`) |
| `Comment` | *(optional)* Shown in the per-client summary for samples below the low-read threshold, and in the QC report |
| `Location` | *(optional, QC report only)* Plate position; carried into the QC report |
| `Sample_Type` | *(optional, QC report only)* `SAMPLE` (default when blank), `NTC`, `BLANK` or `POS_CTRL`. Controls get flag `CTRL` and are left out of run medians |
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
filter_stage=post                  # post | pre -- see "Filter stage" above
run_qc_report=1                    # 1 = write qc_report/ at the end of the job

# --- Variables index, change as needed ---
min_quality=15                     # Chopper -q
min_length=400                     # Chopper --minlength (ITS default)
max_length=1200                    # Chopper --maxlength (ITS default)

primer_mismatch=9                  # Minibar -E
barcode_mismatch=1                 # Minibar -e
```

**Both** barcode references are set explicitly and passed to `generate_primer_setup.py` with `-b`, so nothing depends on the generator's built-in default. Both files must cover every barcode ID used in the sample sheet's `Barcode` column.

If you swap in a different barcode set, also re-check `min_length` / `max_length` and `primer_mismatch` -- both are tuned for the IDT ITS construct. `qc_config.yaml` mirrors `min_length` / `max_length` (`filter.minknow_min_len` / `filter.chopper_max_len`); change both together, or the QC report warns that the windows differ.

Everything in the `DONT-CHANGE` section below resolves automatically from the variables above: `raw_input_dir` is found from `run_name`, and the two generated primer setup files are derived from `primer_output_dir` + `date_stamp`.

The wrapper points Minibar's per-file working directories at `$PBS_JOBFS` (node-local disk) and passes `--minibar-work-dir` to the core script; `pre` mode's Chopper output goes there too. Only the merged per-sample FASTQs are written to `output_dir`. Nothing outside `/g/data/vz35` and `$PBS_JOBFS` is touched (no `/tmp`, no `/scratch`), so `-l storage=gdata/vz35` remains correct in both modes.

---

## Usage

### 1. Prepare the sample sheet

Fill in the CSV with `Client`, `Sample_ID` and `Barcode` columns.
Include a control sample with a known expected yield to validate the run, and mark controls in `Sample_Type` so the QC report treats them as controls.

Optional `Comment` can be put based on the amplification tendency.

### 2. Update the Variables section

Open `run_script.qsub` and set `run_name`, `samplesheet`, `output_dir` and `filter_stage`.
If the barcode set has changed, point `barcode_file_minibar` and `barcode_file_cutadapt` at the new pair and re-check `min_length` / `max_length` and `primer_mismatch`.

### 3. Submit the job

```bash
qsub run_script.qsub
```

PBS resources requested: 16 CPUs, 64 GB RAM, 200 GB jobfs, 8 h walltime, project `vz35`, queue `biodev`, storage `gdata/vz35`.

These are the 16S workflow's v6.00 values, which uses the same parallel Minibar step; they have **not yet been confirmed on an ITS run**:

- **ncpus=16 / mem=64GB.** Step 2 runs Minibar on several `fastq_pass` chunks at once (one subprocess per worker, `$PBS_NCPUS` workers). The speed-up is bounded by the number of large chunks and then by the single largest one. Check the PBS resource footer for peak memory on the first run and trim if it is far below 64 GB.
- **jobfs is 200 GB** (was 10 GB). In `post` mode Minibar's per-file working directories hold unfiltered, uncompressed reads and are much larger than in `pre` mode. On queues that ration jobfs per core, 200 GB may need more than a few CPUs to be schedulable — raise `ncpus` before lowering `jobfs`.
- **Walltime is 8 h** (was 20 h for the serial loop). `post` mode feeds Minibar every read in `fastq_pass`, including the > 1.2 kb ones, and Minibar is the dominant cost, so check the step 2 timestamps in the job log on the first `post` run before trusting 8 h.

### 4. Validate with the control run

After the job completes, open `qc_report/qc_report.html`, `run_log_<date>/read_counts_summary.txt` and the control sample's client `summary.txt` inside `integrated_demultiplexing/<Client>/`.
Confirm that the control sample's read count and percentage match the expected values for that barcode, and check `run_log_<date>/cutadapt_summary.txt` for a plausible kept / discarded / reverse-complemented split.
Once validated, the same script can be resubmitted for any other sample set by updating the Variables section.

If a sample came out low, `read_counts_summary.txt` now says *where* the reads went — `pct_retained` for filter loss, `pct_pass1_kept` for non-specific amplification — and `demux_qc/demux_qc_summary.tsv` gives the length and quality profile behind those numbers.

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
| `ITS_primer_setup_<date>.txt` | `barcode_file_cutadapt` | Step 5 (Cutadapt) |

Both are written to `sample_sheet/` and passed to the core script as `--minibar-primer-file` and `--primer-file` respectively. If `--minibar-primer-file` is omitted, Minibar falls back to `--primer-file`.

---

### Step 1: Chopper -- quality and length filtering (`pre` mode only)

Skipped entirely when `filter_stage=post`, where raw `fastq_pass/*.fastq.gz` go straight to Minibar and filtering happens at step 4 instead.

Each `*.fastq.gz` file in `fastq_pass/` is piped through Chopper, and reads that pass all three filters are retained:

| Filter | Parameter | Value |
|--------|-----------|-------|
| Quality score | `-q` | >= 15 |
| Minimum length | `--minlength` | 400 bp |
| Maximum length | `--maxlength` | 1,200 bp |

All three are set in the **Variables** section of `run_script.qsub` (`min_quality`, `min_length`, `max_length`) and passed through to the core script; the same values are used at step 4 in `post` mode.

The length window must bracket the expected amplicon for the primer pair in use, so re-check it whenever the primers change. The 400-1,200 bp default is sized for ITS: the main peak sits at ~0.78-0.80 kb with a shoulder at ~0.60-0.70 kb, and concatemers/junk appear above ~1.5 kb.

The equivalent command run per file:

```
zcat <file>.fastq.gz | chopper -q 15 --minlength 400 --maxlength 1200 > <file>_filtered.fastq
```

Filtered files are written **uncompressed** to `$PBS_JOBFS/chopper_filtered/` (falling back to `output_dir/chopper_filtered/` when `PBS_JOBFS` is unset) and deleted automatically after demultiplexing. They are intermediate — Minibar reads them once — and gzipping them at Python's default level 9 cost more than the rest of the step.

---

### Step 2: Minibar -- demultiplexing

Minibar demultiplexes each input file into per-sample FASTQs using the **Minibar** primer setup file.
Each file is processed in its own subdirectory (on `$PBS_JOBFS`), several files at a time (`$PBS_NCPUS` workers), and per-sample results are merged across all `fastq_pass` files at the end in input order, so the merged FASTQs are identical to a serial run. Minibar's own output for each file goes to `minibar_logs/<file>.log` (moved into `run_log_<date>/` at step 7). Zero-read input files are skipped rather than crashing Minibar.

Flags are identical in both filter stages; only the input differs (Chopper-filtered files in `pre`, raw `fastq_pass` files in `post`).

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
Per-file Minibar subdirectories are removed after merging.

---

### Step 3: Per-sample demux QC

`demux_qc.py` makes a single streaming pass over each merged `sample_<SampleID>.fastq` **before** any filtering, and reports:

| Metric | Detail |
|--------|--------|
| Raw demuxed read count | Reads Minibar assigned to that barcode |
| Length histogram | 50 bp bins, plus exact N50, median, and modal length |
| Mean read Q | Per-read mean error probability converted to Phred, then averaged over reads |
| Q histogram | Integer bins of that per-read mean Q, for all reads and for reads inside `[min_length, max_length]`; `median_q` = median bin of the in-window reads. For display only -- read quality is not flagged |
| % within window | Share of reads inside `[min_length, max_length]` |

`sample_unk` and `sample_Multiple_Matches` are profiled too — an inflated unassigned bin is the first sign of a demux problem, and its length profile says whether the cause is chemistry or barcode matching.

Outputs, written to `demux_qc/`:

| File | Contents |
|------|----------|
| `<SampleID>.tsv` | One row per 50 bp bin; scalar stats in a `#`-prefixed header block (read with `pandas.read_csv(..., sep='\t', comment='#')`) |
| `<SampleID>.q.tsv` | One row per integer mean-Q bin: `q_bin`, `count_all`, `count_in_window`; same `#` header convention |
| `demux_qc_summary.tsv` | One row per sample, all scalar stats (incl. `median_q`), real samples first |

No read is altered.

The step runs in **both** filter stages, but means different things in each: in `post` mode it profiles genuinely raw demultiplexed reads; in `pre` mode those reads have already been through Chopper. Every QC output records which stage produced it, so the two are never confused. Skip with `--no-demux-qc`.

---

### Step 4: Chopper -- per-sample filtering (`post` mode only)

Skipped when `filter_stage=pre`, where Chopper already ran at step 1.

Applies the same `-q` / `--minlength` / `--maxlength` filter as step 1, but to each merged `sample_<SampleID>.fastq` individually, filtering in place. Pre-filter counts are carried over from step 3 so the demux loss and the filter loss stay separable in the run summary.

`sample_unk` and `sample_Multiple_Matches` are deliberately **left unfiltered** — their raw length profile is the whole reason they are diagnostic. They are therefore also excluded from the TOTAL `pct_retained` denominator in `read_counts_summary.txt`.

---

### Step 5: Cutadapt -- orientation normalisation and primer trimming

Run on each `sample_<SampleID>.fastq` after merging (and, in `post` mode, after step 4), using the **cutadapt** primer setup file (full construct). Reads from `sample_unk` and `sample_Multiple_Matches` are left untouched, as are samples with no matching row in the primer setup file.

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

### Step 6: Organise by client and summarise

Demultiplexed files are moved into `integrated_demultiplexing/<Client>/` based on the sample sheet. If cutadapt skipped a sample, the pre-cutadapt `sample_<id>.fastq` is moved instead and flagged in the log as `Moved (uncut)`.

| File | Contents |
|------|----------|
| `integrated_demultiplexing/<Client>/summary.txt` | Per-sample read counts for that client; shows the `Comment` column for samples below 15,000 reads |
| `integrated_demultiplexing/<Client>/<Client>_qc_report.html` | Client-safe per-sample QC report, written by `qc_report.py` (see "Per-sample QC report") |
| `read_counts_summary.txt` | Run-level totals plus the per-stage read-count breakdown below |

`read_counts_summary.txt` tracks the read count through every stage, so the two-stage loss is separable. Its header records ` Filter stage: post|pre`.

| Column | Meaning |
|--------|---------|
| `Reads` | Final read count in the delivered FASTQ |
| `% of Input` | Share of total raw (`post`) or total filtered (`pre`) input |
| `raw_demuxed` | Reads Minibar assigned, before filtering |
| `post_chopper` | Reads surviving the Chopper length/Q filter |
| `pct_retained` | `post_chopper / raw_demuxed` |
| `post_cutadapt_pass1` | Reads with a findable forward construct |
| `pct_pass1_kept` | `post_cutadapt_pass1 / cutadapt input` |

**`pct_pass1_kept` is the non-specificity signal.** A barcode with a healthy `raw_demuxed` but a low `pct_pass1_kept` means Minibar found the right barcode on reads that are not the expected product — correct barcode, wrong amplicon. A `post_chopper` count alone hides this.

In `pre` mode no raw demuxed count exists (reads were already filtered when Minibar saw them), so `raw_demuxed` and `pct_retained` are left **blank** rather than being filled with a post-filter number.

`SUCCESSFULLY DEMULTIPLEXED` is computed as total input minus `sample_unk` minus `sample_Multiple_Matches`.

The low-read threshold is the `LOW_READ_THRESHOLD` constant in `filter_chopper_demux_minibar.py` (set to `None` to disable).

### Step 7: log collection

All top-level `*.log` and `*.txt` files (including `read_counts_summary.txt` and `cutadapt_summary.txt`) and the `cutadapt_logs/` and `minibar_logs/` directories are moved into `run_log_<date>/` to keep the output root clean. Per-client `summary.txt` files stay with their client directory. `demux_qc/` is a directory of `.tsv` files and stays where it is.

---

## Output structure

```
minibar_output/ONT_ITS_TBC_<date>/
├── integrated_demultiplexing/
│   ├── <ClientA>/
│   │   ├── <SampleID>.fastq.gz
│   │   ├── summary.txt
│   │   └── <ClientA>_qc_report.html
│   ├── <ClientB>/
│   │   └── ...
│   ├── sample_Multiple_Matches.fastq
│   └── sample_unk.fastq
├── demux_qc/
│   ├── <SampleID>.tsv
│   ├── <SampleID>.q.tsv
│   ├── unk.tsv
│   ├── Multiple_Matches.tsv
│   └── demux_qc_summary.tsv
├── qc_report/                       # see "Per-sample QC report"
└── run_log_<date>/
    ├── read_counts_summary.txt
    ├── cutadapt_summary.txt
    ├── cutadapt_logs/
    │   ├── <SampleID>.pass1_orient.log
    │   └── <SampleID>.pass2_trim3.log
    └── minibar_logs/
        └── <fastq_pass file>.log
```

- `integrated_demultiplexing/` -- reads merged across all `fastq_pass` files and split by client
- `sample_unk.fastq` -- reads that did not match any barcode
- `sample_Multiple_Matches.fastq` -- reads that matched more than one barcode
- `summary.txt` -- per-client read counts; low-count samples (< 15,000 reads) show the `Comment` column from the sample sheet
- `demux_qc/` -- per-sample length histograms and quality stats taken before filtering
- `read_counts_summary.txt` -- overall run summary; ends up in `run_log_<date>/` after step 7

---

## Per-sample QC report

`qc_report.py` produces one table per run with one row per sample, a PASS / CHECK / FAIL flag and reason codes, plus run-level verdicts and length/Q heatmaps. It judges only what the facility is responsible for -- did the sample amplify and sequence into usable reads. Chimera %, classified / unclassified reads and NTC read ratios are the client's responsibility and are not reported.

**Thresholds are provisional.** `qc_config.yaml` (`config_version: 1.0-its`) starts from the 16S values and the expected ITS length profile; none has been calibrated on ITS data yet. Calibrate on past ITS runs including **fungal mock / positive controls** — not ZymoBIOMICS D6305, which is a bacterial 16S standard. Until then treat every flag as a prompt to look, not a verdict.

**Read quality is not flagged.** Every delivered read already passed MinKNOW Q15. The Q distribution is still shown (heatmap and drill-down), and `median_q` is reported, for information only.

The main output is **`qc_report/qc_report.html`**: a single self-contained page (no network access needed -- open it from Gadi, a shared drive or an email attachment) with the run summary, interactive length and Q heatmaps, a click-to-drill-down view per sample and a sortable table. `integrated_demultiplexing/<Client>/<Client>_qc_report.html` is the same page cut to one client; it sits in the client's delivery folder so it ships with the reads.

**Client-safe HTML.** No HTML report contains internal paths or identifiers (no `/g/data`, `/scratch`, `/home`, `vz35`, usernames, sample-sheet or config paths); the config is identified by version and short hash, and paths in warnings are cut to file names. A per-client report contains nothing about other clients -- no sample IDs, barcodes, comments or client names, not even in the embedded data. Its header counts are recomputed for that client; run-level flags appear only as anonymous text (e.g. `RUN: LONG (median 18%, FAIL)`), the drill-down's run-median line is an anonymous aggregate, and warnings and catch-all bin counts are left out.

It is a pure reader of `demux_qc/` and `read_counts_summary.txt`: no FASTQ access, no Minibar, no Chopper. It takes seconds and runs automatically at the end of `run_script.qsub` (`run_qc_report=1`). A QC failure there is reported but does not fail the job. To re-run it -- after changing thresholds, on a past run, or if the in-job report failed:

```bash
qsub -v output_dir=<run output dir>,samplesheet=<sheet.csv> run_qc_report.qsub
# or by hand
python3 qc_report.py --output-dir <run output dir> --samplesheet <sheet.csv>
```

Pass `output_dir` explicitly: `run_script.qsub` stamps it with the date its job *starts*, which is not necessarily the day you submit.

### Metrics (all exact; four flagged for ITS)

| Metric | Source | Better | Catches | ITS |
|---|---|---|---|---|
| `delivered_reads` | `Reads` in `read_counts_summary.txt` | higher | Depth vs the 15,000 floor | flagged |
| `primer_found` | `pct_pass1_kept` | higher | Correct barcode, wrong amplicon | flagged |
| `long_frac` | 1 − `pct_retained` (post mode only) | lower | Chimeras, concatemers, ligation artefacts (> 1.2 kb) | flagged |
| `core_frac` | `demux_qc` length histogram: reads in `core_window` (450–1,000 bp) / reads in 400–1,200 bp | higher | Partial / off-target products outside the expected ITS range | flagged |
| `modal_dev_bp` | \|modal length − run median modal length\| | lower | Wrong product / assay mix-up | **info only** |

Informational, never flagged: `median_q`, the median mean-read-Q bin from the `demux_qc` Q histogram, and — for ITS — `modal_dev_bp`.

**Modal length is informational for ITS.** ITS length legitimately varies with community composition (ITS1/ITS2 length differs between fungal taxa), so a sample dominated by short- or long-ITS taxa sits far from the run's median modal length without anything being wrong. The ~650 bp shoulder can even overtake the ~790 bp peak and make the modal length jump by ~140 bp between otherwise similar samples. `modal_dev_bp` and `modal_length` are therefore reported (table, drill-down, `run_qc.txt`) but carry `info_only: true` in `qc_config.yaml` and never set a flag, a reason code or a run flag.

`long_frac = 1 − pct_retained` measures the > 1.2 kb fraction only because MinKNOW's 400 bp floor on ITS runs equals Chopper's `--minlength`, so Chopper removes nothing at the short end. If either changes, `long_frac` would also count short-read loss.

`core_frac` and `median_q` are computed on reads within 400–1,200 bp, so they describe the delivered population. `raw_demuxed` in the table is Minibar's count before Chopper.

**Known limit (ITS blind spot).** Anything shorter than `minknow_min_len` (400 bp) is discarded by MinKNOW and never reaches the pipeline. Primer dimers and short non-specific products therefore cannot be measured per sample; a short-product failure can only show up as low `delivered_reads` or a low `core_frac`. Short and low-Q fractions are not reported.

**`core_window` is provisional.** 450–1,000 bp is wide enough to keep both the ~790 bp peak and the ~650 bp shoulder. Confirm it against the `demux_qc/<SampleID>.tsv` length histograms (or the length heatmap) of a real ITS run before relying on `core_frac`.

**Pre-mode runs.** Chopper ran before Minibar, so no > 1.2 kb reads were demultiplexed: `long_frac` is `NA (pre mode)` for every sample and `run_qc.txt` says so at the top. The other metrics still work. Runs whose `demux_qc/` predates the Q histogram get `median_q = NA`, with a warning.

### Flagging

Per metric 0 = pass, 1 = check, 2 = fail; the sample's flag is the worst flagged metric.

- Absolute thresholds from `qc_config.yaml` by default.
- If a metric's **run median** is already at CHECK or worse, a run flag is raised (e.g. `RUN: LONG (median 18%, FAIL)`) and samples are judged on that metric **relative to the run**: robust z = (value − median) / max(1.4826 × MAD, `mad_floor`), CHECK at z ≥ 3, FAIL at z ≥ 5, in the bad direction only. One bad flowcell therefore raises one run flag instead of turning the whole table red. Fewer than `min_samples` real samples keeps absolute thresholds.
- `delivered_reads` is `absolute_only`: 15,000 reads is a deliverable floor.
- `info_only` metrics (`modal_dev_bp` for ITS) get no level at all: they are never part of the flag, the reasons or a run flag, and their run median is listed as `info only` in `run_qc.txt`.
- Samples with fewer than `min_reads_for_stats` (200) raw demultiplexed reads get `NA (n<200)` in every metric except `delivered_reads`, which is still flagged.
- Controls (`Sample_Type` NTC / BLANK / POS_CTRL) are reported with flag `CTRL` and left out of every run median. A control with an empty `Client` or `Client = ALL` is shared: it appears in every per-client report and TSV. A control with a named client goes only to that client. (The pipeline itself skips rows with an empty `Client`, so such a control has no delivered reads; use `ALL` if it should be delivered.)
- Flags are always run-based: run medians, relative-mode thresholds and the modal-length median come from all real samples of the run. Per-client reports only filter rows, so a sample has the same flag in every report.

Reason codes list every non-pass flagged metric, worst first, e.g. `FAIL: LOW_CORE(48%);LONG(22%)`; relative-mode entries carry their z (`PRIMER(30%,z=17.8)`). Codes: `READS`, `PRIMER`, `LONG`, `LOW_CORE` (and `MODAL`, only if `modal_dev_bp` is made a flagged metric again).

**PASS / CHECK / FAIL is triage, not sign-off.** Review `qc_report.html` (or `run_qc.txt` and the heatmaps) before release. Bump `config_version` whenever a threshold changes -- it is recorded in every output.

### `qc_config.yaml`

Everything assay-specific lives here, so the same `qc_report.py` serves 16S and ITS:

| Key | ITS value | Meaning |
|---|---|---|
| `assay` / `title` | `ITS` / `ITS QC report` | Report title and headers (default `16S` / `16S QC report` when absent) |
| `filter.minknow_min_len` | 400 | MinKNOW minimum read length on ITS runs; lower edge of the delivered window |
| `filter.minknow_min_q` | 15 | MinKNOW pass Q |
| `filter.chopper_max_len` | 1200 | Chopper `--maxlength`; upper edge of the delivered window, and what "long" means (`> 1.2 kb` in every label) |
| `core_window` | 450–1,000 | `core_frac` numerator; multiples of 50, inside the delivered window |
| `plot.len_min` / `plot.len_max` | 400 / 2,000 | Length axis of the heatmaps and drill-down; must start at or below `minknow_min_len` and end above `chopper_max_len` (default `minknow_min_len` .. 2 × `chopper_max_len`) |
| `metrics.<m>.info_only` | `true` for `modal_dev_bp` | Report but never flag; `check` / `fail` not needed |

### Outputs -- `<output_dir>/qc_report/`

| File | Contents |
|---|---|
| `qc_report.html` | **Main review view** (internal, all samples). One self-contained file (data embedded as JSON, inline JS/CSS, no external libraries or network). Header with filter stage, config version, flag counts and run-level flags; filter by client / flag / sample name and sort by any column; length and Q heatmaps side by side (hover for details); click a row for that sample's length and Q histograms against the run median, each metric's PASS / CHECK / FAIL (or *info*) badge and threshold, reasons and Comment; sortable table of all samples. Written with the standard library only, so it is produced even without matplotlib |
| `per_sample_qc.tsv` | One row per sample: client, sample, barcode, type, `raw_demuxed`, every metric, `modal_length`, `median_q` (info), flag, reasons, Comment. `#` header records assay, config version, filter stage, windows and which metrics are info-only |
| `per_client/<Client>_qc.tsv` | Per-client cut of the same table, plus shared controls |
| `../integrated_demultiplexing/<Client>/<Client>_qc_report.html` | Per-client cut of `qc_report.html` (same layout and interactivity), written into the client's delivery folder next to its reads and `summary.txt`; see "Client-safe HTML" above. If the run dir has no `integrated_demultiplexing/`, it goes to `qc_report/per_client/` with a warning |
| `run_qc.txt` | Assay, config version, filter stage, windows, run-level flags and medians, flag counts, `unk` / `Multiple_Matches` share, warnings |
| `qc_bins.tsv.gz` | Long format (sample, kind, bin, count); kind `len` (50 bp, all reads), `q` (in-window reads), `q_all` -- for re-plotting / re-thresholding |
| `length_heatmap.png`, `q_heatmap.png` | Static copies for email / records. Rows = samples (FAIL, CHECK, PASS, CTRL, then sheet order), row-normalised. Length axis 400–2,000 bp shows the > 1.2 kb tail, core window shaded, dashed lines at 400 and 1,200 bp; the Q heatmap is informational and has no threshold line. Flag strip on the left, flagged rows labelled |
| `flagged/<Sample>.png` | Only with `--png-drilldown`: length and Q histograms for CHECK / FAIL samples (the HTML report has the same drill-down) |

### Requirements and tests

Needs PyYAML, and matplotlib for the static PNGs (the tables and `qc_report.html` are written without it). Both come from the shared venv `/g/data/vz35/ONT_16s_workflow/tools/qc_env`; `run_script.qsub` and `run_qc_report.qsub` call its python directly (`qc_python=.../qc_env/bin/python3`) after `module load python3/3.12.1`, so the core pipeline's environment is untouched. The wrapper checks for both packages before starting. To (re)create the venv:

```bash
module load python3/3.12.1
python3 -m venv /g/data/vz35/ONT_16s_workflow/tools/qc_env
/g/data/vz35/ONT_16s_workflow/tools/qc_env/bin/pip install pyyaml matplotlib
```

`tools_dir` (and `qc_python`, `qc_config`, `samplesheet`) can be overridden with `qsub -v`, so a dev copy can be tested without editing the scripts: `qsub -v output_dir=<run>,samplesheet=<sheet.csv>,tools_dir=<dev checkout> run_qc_report.qsub`.

`python tests/test_qc_report.py` builds fake ITS runs in post and pre mode (real `demux_qc` + production `read_counts_summary.txt` writer) and checks: the ITS config loads; a core window below 1 kb is accepted and windows outside the delivered range are rejected; a legacy 16S config still loads with 16S defaults; the Q bins and `median_q` against a direct calculation; every metric and flag (Q and `modal_dev_bp` never flagged, even with a wildly spread modal length); the counts parser incl. pre-mode NA handling; relative mode; controls; old-config warnings; that the HTML title, long-read label and length axis come from the config (no `16S` / `2 kb` left); that `qc_report.html` parses, contains every sample, has no `http://` / `https://` references and is still written without matplotlib; and that no HTML report — internal or per-client — contains `/g/data`, `/scratch`, `/home` or `vz35`.

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
| `--filter-stage {post,pre}` | `post` | When Chopper runs; see "Filter stage" |
| `--diagnostic` | off | Shorthand for `--filter-stage post --min-length 100 --max-length 5000`. Explicitly given flags still win |
| `--minibar-work-dir` | `$PBS_JOBFS/minibar_work`, else the output dir | Where Minibar's per-file working directories are created |
| `--no-demux-qc` | off | Skip step 3 |
| `--chopper` | Gadi 16S tools path | Chopper binary |
| `--minibar` | Gadi 16S tools path | `minibar.py` |
| `--cutadapt` | Gadi 16S tools path | Cutadapt binary |
| `--min-quality` | 15 | Chopper `-q` |
| `--min-length` (`--minlength`) | 400 | Chopper `--minlength` |
| `--max-length` (`--maxlength`) | 1200 | Chopper `--maxlength` |
| `--barcode-mismatch` | 1 | Minibar `-e` |
| `--primer-mismatch` | 5 | Minibar `-E` (wrapper passes 9) |
| `--cutadapt-threads` | 4 | Cutadapt `-j` |
| `--cutadapt-error-rate` | 0.2 | Cutadapt `-e` |

### Diagnostic mode

For a first look at an unfamiliar run, `--diagnostic` forces `post` mode with a 100–5,000 bp window so step 3 observes the **real** length distribution instead of only what survives an assumed 400–1,200 bp window. Read `demux_qc/demux_qc_summary.tsv`, pick a window that matches the observed modal length and N50, then do the production run with those values. (MinKNOW's 400 bp floor still applies, so nothing below 400 bp is visible even here.)

### Running demux QC standalone

`demux_qc.py` also works on an existing run directory without rerunning the pipeline:

```bash
python3 demux_qc.py \
    --integrated-dir /path/to/run_output/integrated_demultiplexing \
    [--qc-dir /path/to/demux_qc] \
    [--minlength 400] [--maxlength 1200] [--filter-stage pre|post]
```

It only reads the FASTQs — it never rewrites reads — so it is safe to run against a completed run. Note that after step 6 the per-sample files have been renamed to `<SampleID>.fastq.gz` and moved into client subdirectories, so a standalone run against the top level will only find the `sample_unk` / `sample_Multiple_Matches` bins.
