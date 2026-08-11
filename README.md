# MetaBarFlow

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Made with DADA2](https://img.shields.io/badge/ASVs-DADA2-1f77b4.svg)](https://benjjneb.github.io/dada2/)

A generalised, marker-agnostic metabarcoding pipeline for Illumina paired-end amplicon data.  
Supports 18S (protists), ITS (fungi), 16S (bacteria/archaea), and any other amplicon marker.

---

## Pipeline overview

The workflow has two stages: the **bioinformatics pipeline** (steps 1–5, run via `scripts/main.sh`)
and the downstream **statistics stage** (run in R, see [Statistics stage](#statistics-stage)).

```
Bioinformatics pipeline (scripts/main.sh)
  Step 1 — Merge FastQ        : pool R1/R2 files from multiple sequencing runs per sample
  Step 2 — Quality check      : FastQC + MultiQC on merged reads
  Step 3 — Primer trimming    : remove primers and adapters with Cutadapt
  Step 4 — ASV generation     : filter, denoise, merge, and chimera-filter with DADA2
  Step 5 — Taxonomy assignment: assign taxonomy to ASVs with VSEARCH against a reference DB

Statistics stage (Stats/)
  Import      : load the per-marker MR/ASS tables
  Filtering   : taxonomic / rare-ASV / quality / sample filtering, rarefaction, depth correction
  Diagnostics : dataset-quality figures and tables
```

Each pipeline step writes its output to a subdirectory of `results/<MARKER>/` and appends a row to
`results/<MARKER>/logs/read_tracking.tsv` so you can track read loss at every stage.

---

## Directory structure

```
MetaBarFlow/
├── scripts/
│   ├── main.sh                 ← master launcher — run this
│   ├── 01_merge_fastq.sh
│   ├── 02_quality_check.sh
│   ├── 03_trim_primers.sh
│   ├── 04_DADA2_ASV.R
│   ├── 05_taxonomy.R
│   └── utils.sh
│
├── config/
│   ├── config_template.sh      ← copy + edit to create a new marker config
│   └── config_template.R
│
├── Pipeline_18S/               ← 18S-specific configuration (ready to use)
│   ├── config.sh
│   └── config.R
├── Pipeline_ITS/               ← ITS-specific configuration (ready to use)
│   ├── config.sh
│   └── config.R
├── Pipeline_16S/               ← 16S-specific configuration (ready to use)
│   ├── config.sh
│   └── config.R
│
├── RAW_dataset/                ← input FASTQ files (one subdirectory per run)
│   ├── 18S/
│   │   ├── Heger18SV4_1/
│   │   └── Heger18SV4_2/
│   ├── ITS/
│   │   ├── Heger_ITS_1/
│   │   └── Heger_ITS_2/
│   └── 16S/
│       ├── HegerV4V5_1/
│       └── HegerV4V5_2/
│
├── databases/                  ← reference databases (not included in repo)
│   ├── PR2/                    ← PR2 v5.1.0 (18S)
│   ├── UNITE/                  ← UNITE v9 (ITS)
│   └── SILVA/                  ← SILVA 138.2 (16S)
│
├── results/<MARKER>/          ← created automatically on first run
│   ├── 01_merged/
│   ├── 02_QC/
│   ├── 03_trimmed/
│   ├── 04_ASV/
│   ├── 05_taxonomy/
│   └── logs/
│       └── read_tracking.tsv
│
└── Stats/                      ← statistics stage (run in R, after the pipeline)
    ├── 00_import.R             ← import MR/ASS → DS object
    ├── 01_filtering.R          ← filtering, rarefaction, depth correction
    ├── 02_dataset_quality.R    ← diagnostic figures + tables
    └── README.md
```

> **Note:** Marker configuration lives in `Pipeline_<MARKER>/` directories (one per marker),
> while pipeline output is written to `results/<MARKER>/`. Keep these two separate.

---

## Quick start

### 1. Install dependencies

See [INSTALL.md](INSTALL.md) for full instructions.

### 2. Download reference databases

Place databases in the `databases/` subdirectory:

| Marker | Database | Directory |
|--------|----------|-----------|
| 18S    | PR2 v5.1.0 | `databases/PR2/` |
| ITS    | UNITE v9 (19.02.2025) | `databases/UNITE/` |
| 16S    | SILVA 138.2 NR99 | `databases/SILVA/` |

Download links are listed in INSTALL.md.

### 3. Place raw FASTQ files

Paired-end FASTQ files (`*_R1.fastq.gz` / `*_R2.fastq.gz`) must be placed in a run
directory under `RAW_dataset/<MARKER>/`. Each run directory holds all files for one
sequencing run. Multiple run directories are merged in step 1.

### 4. Run the pipeline

```bash
# Run all 5 steps for 18S
bash scripts/main.sh Pipeline_18S/config.sh --steps 1,2,3,4,5

# Run only steps 1–4 (skip taxonomy — database not present)
bash scripts/main.sh Pipeline_18S/config.sh --steps 1,2,3,4

# Re-run step 3 alone (after adjusting cutadapt parameters)
bash scripts/main.sh Pipeline_18S/config.sh --steps 3

# Run ITS pipeline
bash scripts/main.sh Pipeline_ITS/config.sh --steps 1,2,3,4,5

# Run 16S pipeline
bash scripts/main.sh Pipeline_16S/config.sh --steps 1,2,3,4,5
```

> **Always run from the project root** (`MetaBarFlow/`), not from inside a subdirectory.

### 5. Check outputs

| File | Description |
|------|-------------|
| `results/<MARKER>/logs/read_tracking.tsv` | Read counts at each step |
| `results/<MARKER>/02_QC/multiqc_report.html` | MultiQC quality report |
| `results/<MARKER>/04_ASV/<MARKER>_MR.csv` | ASV abundance matrix |
| `results/<MARKER>/04_ASV/<MARKER>_Fasta.fasta` | ASV sequences |
| `results/<MARKER>/05_taxonomy/<MARKER>_ASS.csv` | ASV table with taxonomy |

---

## Try it on the example data

A small **real** example dataset is bundled under `RAW_dataset/` for all three markers
(4 samples × 2 runs each) so you can verify a fresh installation without supplying your own
data. The files are real sequencing reads **subsampled to ~3,000 read pairs each** (see
[`test/subsample_raw.sh`](test/subsample_raw.sh)) to keep the repository small. 18S is the most
complete and the recommended smoke test:

```bash
conda activate metabarflow
# Steps 1–4 need no reference database; add step 5 once PR2 is in databases/
bash scripts/main.sh Pipeline_18S/config.sh --steps 1,2,3,4
```

Outputs appear under `results/18S/`. See [`test/README.md`](test/README.md) for details,
expected outputs, and the two intentional real-world edge cases (a missing 16S R2 file and
mixed ITS sample naming).

---

## Adding a new marker

```bash
# 1. Copy the template configs
cp config/config_template.sh Pipeline_MYMARKER/config.sh
cp config/config_template.R  Pipeline_MYMARKER/config.R

# 2. Edit both files — replace all PLACEHOLDER values
#    Minimum required: MARKER, primers, RUN_DIRS, DB paths, truncation lengths

# 3. Create the raw data directory and add FASTQ files
mkdir -p RAW_dataset/MYMARKER/Run_1
# copy or symlink FASTQ files into RAW_dataset/MYMARKER/Run_1/

# 4. Run the pipeline
bash scripts/main.sh Pipeline_MYMARKER/config.sh --steps 1,2,3,4,5
```

---

## Configuration reference

All pipeline behaviour is controlled by two marker-specific config files:

- `Pipeline_<MARKER>/config.sh` — Bash settings (primers, cutadapt, paths, threads)
- `Pipeline_<MARKER>/config.R`  — R/DADA2 settings (filter params, taxonomy levels)

Documented templates with explanations for every parameter are in `config/`.

---

## Read tracking

`results/<MARKER>/logs/read_tracking.tsv` records per-sample read counts after each step:

| Column | Description |
|--------|-------------|
| `sample` | Sample name |
| `step` | Pipeline step (1–5) |
| `reads_in` | Input reads |
| `reads_out` | Reads passing this step |
| `pct_retained` | Percentage retained |

---

## Statistics stage

The second stage of the MetaBarFlow workflow. Once the pipeline has produced the per-marker
`MR`/`ASS` tables, the **statistics stage** under [`Stats/`](Stats/) turns them into filtered
community matrices and dataset-quality diagnostic figures (read-depth distribution,
reads-vs-richness, rarefaction curves, ASV length, % identity, filtering journey). It is a
standalone, marker-agnostic, **metadata-free** R workflow, run after the pipeline:

```bash
cd Stats
Rscript 00_import.R && Rscript 01_filtering.R && Rscript 02_dataset_quality.R
```

Requires R packages `vegan`, `ggplot2`, `patchwork`, `scales`. See [`Stats/README.md`](Stats/README.md)
for configuration and outputs.

---

## Known issues

- **16S test data**: `SF_035_S130_L001_R2_001.fastq.gz` is missing from `RAW_dataset/16S/HegerV4V5_1/`. Step 1 will skip this sample with a warning.
- **Step 5 taxonomy**: Requires database files in `databases/`. Step 5 will fail if the database for the selected marker is not present.
- **DADA2 learnErrors multi-threading**: `THREADS_LEARN` is kept at 1 in all configs to avoid a known segfault on Linux with DADA2 < 1.38.

---

## How to cite

If you use MetaBarFlow, please cite the pipeline itself. Citation metadata is in
[`CITATION.cff`](CITATION.cff) — on GitHub, use the **"Cite this repository"** button in the
sidebar. Once a release is archived on Zenodo, cite it via its DOI.

Please also cite the tools MetaBarFlow depends on:

- **DADA2**: Callahan et al. (2016) *Nature Methods* 13:581–583
- **Cutadapt**: Martin (2011) *EMBnet.journal* 17:10–12
- **VSEARCH**: Rognes et al. (2016) *PeerJ* 4:e2584
- **PR2**: Guillou et al. (2013) *Nucleic Acids Research* 41:D597–D604
- **UNITE**: Nilsson et al. (2019) *Nature* 23:49–51
- **SILVA**: Quast et al. (2013) *Nucleic Acids Research* 41:D590–D596
- **vegan** (if you use the statistics stage): Oksanen et al. (2024) *vegan: Community Ecology Package*. R package.

---

## Contributing

Contributions are welcome — bug reports, documentation, and new marker configurations.
See [`CONTRIBUTING.md`](CONTRIBUTING.md) and the [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
Changes are recorded in [`CHANGELOG.md`](CHANGELOG.md).

---

## License

Released under the [MIT License](LICENSE). © 2026 David Singer.
