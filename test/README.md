# MetaBarFlow — Test Dataset

This document describes the test data included with MetaBarFlow and explains how to use it
to verify a fresh installation.

The test FASTQ files are **real sequencing data** (not synthetic), kept small enough to run
a complete pipeline test in minutes.  
They live in `RAW_dataset/` at the project root — not in this `test/` directory.

Each file has been **subsampled to at most ~3,000 read pairs** (its first reads) so the whole
example set is only a few MB — light enough to bundle in the repository while still running the
full pipeline end to end. Files that were already smaller are kept as-is. The exact procedure is
reproducible via [`subsample_raw.sh`](subsample_raw.sh):

```bash
bash test/subsample_raw.sh 3000 RAW_dataset   # cap = read pairs per file
```

---

## Test data layout

```
RAW_dataset/
├── 18S/                        ← RECOMMENDED for testing (most complete)
│   ├── Heger18SV4_1/           run 1 — 4 samples (E_SF_033–036), all R1+R2 present
│   │   ├── E_SF_033_S193_L001_R1_001.fastq.gz
│   │   ├── E_SF_033_S193_L001_R2_001.fastq.gz
│   │   ├── E_SF_034_S205_L001_R1_001.fastq.gz
│   │   ├── E_SF_034_S205_L001_R2_001.fastq.gz
│   │   ├── E_SF_035_S217_L001_R1_001.fastq.gz
│   │   ├── E_SF_035_S217_L001_R2_001.fastq.gz
│   │   ├── E_SF_036_S229_L001_R1_001.fastq.gz
│   │   └── E_SF_036_S229_L001_R2_001.fastq.gz
│   └── Heger18SV4_2/           run 2 — 4 samples (E_SF_033–036), all R1+R2 present
│       ├── E_SF_033_S193_L001_R1_001.fastq.gz
│       ├── E_SF_033_S193_L001_R2_001.fastq.gz
│       ├── E_SF_034_S200_L001_R1_001.fastq.gz
│       ├── E_SF_034_S200_L001_R2_001.fastq.gz
│       ├── E_SF_035_S207_L001_R1_001.fastq.gz
│       ├── E_SF_035_S207_L001_R2_001.fastq.gz
│       ├── E_SF_036_S214_L001_R1_001.fastq.gz
│       └── E_SF_036_S214_L001_R2_001.fastq.gz
│
├── ITS/
│   ├── Heger_ITS_1/            run 1 — 4 samples (F_SF_033–036), all R1+R2 present
│   └── Heger_ITS_2/            run 2 — mixed sample names (F-SF-*, SF-*)
│
└── 16S/
    ├── HegerV4V5_1/            run 1 — 4 samples; SF_035 R2 MISSING (known issue)
    │   ├── SF_033_S112_L001_R1_001.fastq.gz
    │   ├── SF_033_S112_L001_R2_001.fastq.gz
    │   ├── SF_034_S121_L001_R1_001.fastq.gz
    │   ├── SF_034_S121_L001_R2_001.fastq.gz
    │   ├── SF_035_S130_L001_R1_001.fastq.gz   ← R1 only; R2 missing
    │   ├── SF_036_S139_L001_R1_001.fastq.gz
    │   └── SF_036_S139_L001_R2_001.fastq.gz
    └── HegerV4V5_2/            run 2 — 3 samples (SF-36_FB, SF-66, SF-67)
```

---

## Sample naming

Files follow the standard Illumina naming convention:

```
<SampleID>_<SampleIndex>_<Lane>_<Read>_<Set>.fastq.gz
e.g.  E_SF_033_S193_L001_R1_001.fastq.gz
      └──────── ──── ──── ── ───
      SampleID  Idx  Lane  R  Set
```

Step 1 (`01_merge_fastq.sh`) derives sample names automatically by stripping the Illumina
suffix pattern from the R1 filename. The resulting sample name for the file above would be
`E_SF_033`.

---

## Recommended test run: 18S (steps 1–4)

The 18S dataset is the most complete and the best choice for verifying the installation.
Both run directories have all 4 samples with paired R1 and R2 files.

Step 5 (taxonomy) requires the PR2 v5.1.0 database in `databases/PR2/` — see INSTALL.md.
If the database is present, add step 5 to the command below.

```bash
# From the MetaBarFlow project root:
conda activate metabarflow
bash scripts/main.sh Pipeline_18S/config.sh --steps 1,2,3,4
```

Expected outputs after a successful run:

```
results/18S/
├── 01_merged/          merged FASTQ files, one per sample
├── 02_QC/              FastQC + multiqc_report.html
├── 03_trimmed/         primer-trimmed FASTQ files
├── 04_ASV/
│   ├── 18S_Fasta.fasta         ASV sequences
│   ├── 18S_MR.csv              ASV abundance matrix (samples × ASVs)
│   └── seqtab_nochim.rds       DADA2 sequence table (R object)
└── logs/
    ├── main.log
    └── read_tracking.tsv       per-sample read counts at each step
```

---

## Known issues with test data

| Marker | Issue | Impact |
|--------|-------|--------|
| 16S | `SF_035_S130_L001_R2_001.fastq.gz` missing from `HegerV4V5_1/` | Step 1 skips sample SF_035 from run 1; pipeline continues with remaining samples |
| ITS | Run 2 contains mixed naming patterns (`F-SF-*` and `SF-*`) | Step 1 assigns the correct sample names; no pipeline failure expected |

---

## Re-running after a failed step

Each step is idempotent — re-running overwrites previous output for that step only.

```bash
# Re-run only step 3 (e.g. after adjusting cutadapt parameters in config.sh)
bash scripts/main.sh Pipeline_18S/config.sh --steps 3

# Re-run steps 3 and 4 together
bash scripts/main.sh Pipeline_18S/config.sh --steps 3,4
```

---

## Next: statistics stage

Once the pipeline has produced `results/<MARKER>/04_ASV/*_MR.csv` and
`results/<MARKER>/05_taxonomy/*_ASS.csv`, run the statistics stage
(import → filtering → dataset-quality diagnostics) on the example data:

```bash
cd Stats
Rscript 00_import.R && Rscript 01_filtering.R && Rscript 02_dataset_quality.R
```

See [`../Stats/README.md`](../Stats/README.md). The bundled example data is shallow, so the stage
ships with low `RARE_THRESHOLDS` defaults tuned for it.
