# MarkerFlow — Installation Guide

## Requirements

| Tool | Version tested | Role |
|------|---------------|------|
| bash | ≥ 4.0 | pipeline runner |
| R | ≥ 4.3.0 | DADA2 steps |
| cutadapt | 4.9 | primer trimming |
| fastqc | 0.12.1 | per-sample QC |
| multiqc | 1.33 | aggregate QC report |
| vsearch | 2.30.4 | taxonomy assignment |
| DADA2 (R package) | 1.38.0 | ASV generation |
| vegan, ggplot2, patchwork, scales (R packages) | — | statistics stage (`Stats/`) |

All tools except R/DADA2 are available via conda/mamba.

---

## Option A — Conda (recommended)

### 1. Install Miniconda or Mambaforge

```bash
# Miniconda (Linux x86_64)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2. Create the MarkerFlow environment

```bash
conda create -n markerflow -c bioconda -c conda-forge \
    python=3.11 \
    cutadapt=4.9 \
    fastqc=0.12.1 \
    multiqc=1.33 \
    vsearch=2.30.4 \
    r-base=4.3.3

conda activate markerflow
```

### 3. Install DADA2 inside R

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")
```

Verify the installation:

```r
library(dada2)
packageVersion("dada2")   # should print 1.38.0
```

### 3b. R packages for the statistics stage

The `Stats/` statistics stage needs a few CRAN packages. Install them via conda:

```bash
conda install -n markerflow -c conda-forge r-vegan r-ggplot2 r-patchwork r-scales
```

or from within R:

```r
install.packages(c("vegan", "ggplot2", "patchwork", "scales"))
```

### 4. Verify all tools

```bash
conda activate markerflow
cutadapt --version     # 4.9
fastqc --version       # FastQC v0.12.1
multiqc --version      # multiqc, version 1.33
vsearch --version      # vsearch v2.30.4
Rscript -e 'library(dada2); cat(as.character(packageVersion("dada2")), "\n")'
```

---

## Option B — Manual installation

Install each tool individually following their official documentation:

- **cutadapt**: `pip install cutadapt==4.9`
- **FastQC**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: `pip install multiqc==1.33`
- **VSEARCH**: https://github.com/torognes/vsearch/releases
- **R**: https://cran.r-project.org/
- **DADA2**: see step 3 above

---

## Reference databases

Databases must be downloaded separately and placed in the `databases/` subdirectory.
They are **not** included in the repository.

### PR2 v5.1.0 (18S)

```bash
mkdir -p databases/PR2
cd databases/PR2

# mothur-formatted FASTA (used by VSEARCH in step 5)
wget https://github.com/pr2database/pr2database/releases/download/v5.1.0/pr2_version_5.1.0_SSU_mothur.fasta.gz
gunzip pr2_version_5.1.0_SSU_mothur.fasta.gz

# DADA2-formatted FASTA (used by assignTaxonomy in step 4/5)
wget https://github.com/pr2database/pr2database/releases/download/v5.1.0/pr2_version_5.1.0_SSU_dada2.fasta.gz
gunzip pr2_version_5.1.0_SSU_dada2.fasta.gz

cd ../..
```

Expected files after download:
```
databases/PR2/pr2_version_5.1.0_SSU_mothur.fasta
databases/PR2/pr2_version_5.1.0_SSU_dada2.fasta
```

### UNITE v9 (ITS, release 19.02.2025)

```bash
mkdir -p databases/UNITE
cd databases/UNITE

# Download the general release dynamic FASTA from https://unite.ut.ee/repository.php
# File: sh_general_release_dynamic_19.02.2025.tgz
# (direct URL requires accepting UNITE terms — download manually from the site)

# After download and extraction, expected files:
# databases/UNITE/sh_general_release_dynamic_19.02.2025.fasta
# databases/UNITE/sh_general_release_dynamic_19.02.2025_dada2_fixed.fasta

cd ../..
```

### SILVA 138.2 NR99 (16S)

```bash
mkdir -p databases/SILVA
cd databases/SILVA

# Full-length reference (used by VSEARCH)
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz

# DADA2 train set
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.2_train_set.fa.gz
# (keep gzipped — DADA2 reads .gz natively)

cd ../..
```

Expected files after download:
```
databases/SILVA/SILVA_138.2_SSURef_NR99_tax_silva.fasta
databases/SILVA/silva_nr99_v138.2_train_set.fa.gz
```

---

## Windows / WSL2 notes

MarkerFlow is developed and tested on WSL2 (Ubuntu). All commands above work unchanged.

- Use a WSL2 path (`/mnt/d/...`) or a native Linux path — avoid Windows-style paths in bash.
- Ensure `conda activate` works from WSL2: add `conda init bash` to your `~/.bashrc`.
- DADA2 multi-threading on WSL2: `THREADS_LEARN` must remain 1 (set in all provided configs).

---

## Quick check after installation

From the project root, verify the pipeline can be sourced without errors:

```bash
conda activate markerflow
source Pipeline_18S/config.sh
echo "MARKER=${MARKER}, PROJECT_ROOT=${PROJECT_ROOT}"
```

Expected output:
```
MARKER=18S, PROJECT_ROOT=/path/to/MarkerFlow
```
