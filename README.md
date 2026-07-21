# ONT ITS Workflow - BRF @ Gadi

Pipeline for filtering, demultiplexing, orientation normalized and read-counting Oxford Nanopore ITS amplicon data on the NCI Gadi HPC.

---

## Overview

The workflow runs in a single PBS job (`run_script.qsub`) that calls two Python scripts in sequence:

| Stage | Script / Tool | What it does |
|-------|--------------|--------------|
| Pre-step | `generate_primer_setup.py` | Converts the sample sheet CSV into a per-sample primer setup file |
| Step 1 | `filter_chopper_demux_minibar.py` -> **Chopper** | Quality- and length-filters every raw `fastq.gz`  |
| Step 2 | -> **Minibar** | Demultiplexes each filtered file into per-sample FASTQs, then merges across all `fastq_pass` files |
| Step 2.5 | -> **Cutadapt** (two-pass) | Normalises read orientation and trims 5'/3' primers |
| Step 3 | -> organise + summarise | Groups reads by client, writes per-client and run-level read-count summaries |
| Final | -> log collection | Moves all `*.log` / `*.txt` files into `run_log_<date>/` |

After running the pipeline on the control sample and confirming that demultiplexing counts and per-sample yields match expectations, the same `run_script.qsub` can be submitted unchanged for production samples.

---

## Repository contents

```
run_script.qsub                  # PBS job wrapper -- edit Variables section here
generate_primer_setup.py         # Pre-step: sample sheet -> primer setup file
filter_chopper_demux_minibar.py  # Core pipeline: Steps 1, 2, 2.5, 3
cutadapt_2pass.py                # Two-pass cutadapt module (imported by core script)
```

---

## Prerequisites

All tools are expected to be present under `/g/data/vz35/ONT_16s_workflow/`:

| Tool | Default path |
|------|-------------|
| Python 3.12 | via `module load python3/3.12.1` |
| Chopper | `/g/data/vz35/ONT_16s_workflow/tools/chopper/chopper-linux-musl` |
| Minibar | `/g/data/vz35/ONT_16s_workflow/tools/minibar/minibar.py` |
| Cutadapt | `/g/data/vz35/ONT_16s_workflow/tools/cutadapt-env/bin/cutadapt` |
| `generate_primer_setup.py` | `/g/data/vz35/ONT_16s_workflow/tools/brf_ont_16s_workflow/` |
| `filter_chopper_demux_minibar.py` | `/g/data/vz35/ONT_16s_workflow/tools/brf_ont_16s_workflow/` |
| Twist 384 barcode reference | `/g/data/vz35/ONT_16s_workflow/tools/Twist_16S_384_barcode.txt` |

---

## Input files

### 1. Sample sheet (CSV)

Required columns:

| Column | Description |
|--------|-------------|
| `Client` | Client name -- used to organise output into subdirectories |
| `Sample_ID` | Sample identifier (alphanumeric, `-`, `_`; other characters are sanitised) |
| `Barcode` | Twist 384 barcode ID matching a row in the barcode reference |
| `Comment` | *(optional)* Shown in the per-client summary for samples below the low-read threshold |

### 2. Raw reads

PromethION `fastq_pass` directory located at:

```
/g/data/vz35/PromethION_data/sequencer_uploads/<run_name>/
```

The script locates the `fastq_pass` folder automatically from the `run_name` variable.

Only one `fastq_pass` folder is allowed per run.

---

## Configuration

Edit the **Variables** section near the top of `run_script.qsub`:

```bash
run_name=ONT_16S_20260422          # Subdirectory under sequencer_uploads/
samplesheet=/g/data/vz35/ONT_16s_workflow/sample_sheet/16s_samplesheet.csv
output_dir=/g/data/vz35/ONT_16s_workflow/minibar_output/ONT_16S_TBC_<date>

min_quality=15                     # Chopper -q
min_length=400                     # Chopper --minlength (ITS default)
max_length=1200                    # Chopper --maxlength (ITS default)
```

Everything in the `DONT-CHANGE` section below resolves automatically from the first three variables.

---

## Usage

### 1. Prepare the sample sheet

Fill in `16s_samplesheet.csv` with `Client`, `Sample_ID`, and `Barcode` columns.
Include a control sample with a known expected yield to validate the run.

Optional `Comment` can be put based on the amplification tendency.

### 2. Update the Variables section

Open `run_script.qsub` and set `run_name`, `samplesheet`, and `output_dir`.

### 3. Submit the job

```bash
qsub run_script.qsub
```

PBS resources requested: 2 CPUs, 10 GB RAM, 10 GB jobfs, 20 h walltime.

### 4. Validate with the control run

After the job completes, open `read_counts_summary.txt` and the control sample's `summary.txt` inside `integrated_demultiplexing/<Client>/`.
Confirm that the control sample's read count and percentage match the expected values for that barcode.
Once validated, the same script can be resubmitted for any other sample set by updating the Variables section.

### 5. (Optional) Generate the primer setup file standalone

The primer file is generated automatically inside the job. To create it independently:

```bash
python3 generate_primer_setup.py <samplesheet.csv> [-o OUTPUT_DIR] [-b BARCODE_FILE] [--date YYYYMMDD]
```

Writes `ITS_primer_setup_<date>.txt` (tab-separated) with columns:
`SampleID`, `FwIndex`, `FwPrimer`, `RvIndex`, `RvPrimer`.

The default barcode reference is the Twist 384 file. Use `-b` to supply an alternative barcode file for runs with external barcodes.

---

## Pipeline steps in detail

### Pre-step: generate_primer_setup.py

Reads the sample sheet and looks up each barcode in the Twist 384 reference to produce a tab-separated primer setup file.
The file is written to the `sample_sheet/` directory and consumed by all downstream steps.

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

Minibar demultiplexes each filtered file into per-sample FASTQs using the primer setup file.
Each file is processed in its own subdirectory, and per-sample results are merged across all `fastq_pass` files at the end.

| Flag | Value | Meaning |
|------|-------|---------|
| `-e 1` | 1 | Allowed mismatches in barcode |
| `-E 5` | 5 | Allowed mismatches in primer |
| `-l 200` | 200 | Search window at each read end (bp) |
| `-M 2` | 2 | Require barcode match on both ends |
| `-F` | -- | Write each sample to its own file |

**`-M 2` must not be relaxed.** The IDT ITS 96 barcode set is *combinatorial*, not unique-dual: 12 forward indices x 8 reverse indices generate the 96 samples, so each `FwIndex` is shared by 8 samples and each `RvIndex` by 12. Only the *pair* identifies a sample. Under `-M 1`, a single-end match would be ambiguous across 8-12 candidates and reads would be silently cross-assigned between samples instead of landing in `sample_unk`.

`-l 200` is the per-end search window. The barcoded construct is ~52 bp (Illumina tail + 10 bp index + primer), so 200 leaves ample headroom for the few bases of noise ONT reads often carry before the adapter begins. Neither flag changed in the move from 16S to ITS.

Output files: `sample_<SampleID>.fastq`, `sample_unk.fastq`, `sample_Multiple_Matches.fastq`.
Per-file Minibar subdirectories are removed after merging to save storage.

---

### Step 2.5: Cutadapt -- orientation normalisation and primer trimming

Run on each `sample_<SampleID>.fastq` after merging. Reads from `sample_unk` and `sample_Multiple_Matches` are left untouched.

**Pass 1 -- orient and trim 5' primer (strict)**

```
cutadapt -g FWD_PRIMER --revcomp --rename={header} --discard-untrimmed -e 0.2
```

- Searches for the forward primer at the 5' end `-g FWD`
- If the primer is found on the minus strand, the read is reverse-complemented `--revcomp` so all output reads face the same direction
- `--rename={header}` keep the name unchanged
- Reads where the forward primer cannot be found are discarded `--discard-untrimmed` — these are likely noise or off-target sequences

**Pass 2 -- trim 3' primer (tolerant)**

```
cutadapt -a REV_PRIMER_RC -e 0.2
```

- Trims the reverse-complemented reverse primer from the 3' end.
- Reads where the 3' construct is not found are **kept** (truncated reads are retained).

Output: `<SampleID>.fastq.gz` replaces `sample_<SampleID>.fastq`.
Per-sample logs are written to `cutadapt_logs/` and a run-level `cutadapt_summary.txt` is produced.

---

### Step 3: Organise by client and summarise

Demultiplexed files are moved into `integrated_demultiplexing/<Client>/` based on the sample sheet.
Two summary files are written:

| File | Contents |
|------|----------|
| `integrated_demultiplexing/<Client>/summary.txt` | Per-sample read counts for that client; flags samples below 15,000 reads |
| `read_counts_summary.txt` | Run-level totals: total filtered input, total demultiplexed, successfully demultiplexed, and per-sample percentages |

### Final: log collection

All `*.log` and `*.txt` files (including `cutadapt_logs/`) are moved into `run_log_<date>/` to keep the output root clean.

---

## Output structure

```
minibar_output/ONT_16S_TBC_<date>/
├── integrated_demultiplexing/
│   ├── <ClientA>/
│   │   ├── <SampleID>.fastq.gz
│   │   └── summary.txt
│   ├── <ClientB>/
│   │   └── ...
│   ├── sample_Multiple_Matches.fastq
│   └── sample_unk.fastq
├── read_counts_summary.txt
└── run_log_<date>/
    ├── cutadapt_summary.txt
    ├── read_counts_summary.txt
    └── cutadapt_logs/
        ├── <SampleID>.pass1_orient.log
        └── <SampleID>.pass2_trim3.log
```

- `integrated_demultiplexing/` -- reads merged across all `fastq_pass` files and split by client
- `sample_unk.fastq` -- reads that did not match any barcode
- `sample_Multiple_Matches.fastq` -- reads that matched more than one barcode
- `summary.txt` -- per-client read counts; low-count samples (< 15,000 reads) show the `Comment` column from the sample sheet
- `read_counts_summary.txt` -- overall run summary with percentages relative to total filtered input

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
| `-e 1` | 1 | Barcode mismatch tolerance |
| `-E 5` | 5 | Primer mismatch tolerance |
| `-l 200` | 200 | End-window search length (bp) |
| `-M 2` | 2 | Match barcodes on both ends -- required, see Step 2 |
| `-F` | -- | Write separate file per sample |

### Cutadapt

| Parameter | Value | Effect |
|-----------|-------|--------|
| `-e 0.2` | 0.2 | Error rate for adapter matching (both passes) |
| `-j 4` | 4 | Threads per sample |
| `--revcomp` | -- | Pass 1: reverse-complement reads where primer found on minus strand |
| `--discard-untrimmed` | -- | Pass 1: drop reads with no forward primer match |
| `-g FWD` | forward primer | Pass 1: 5'-anchored adapter |
| `-a REV_RC` | RC of reverse primer | Pass 2: 3' adapter |
