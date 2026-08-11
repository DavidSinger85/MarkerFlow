#!/usr/bin/env bash
# =============================================================================
# File      : config_template.sh
# Author    : David Singer
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Description: Template for marker-specific shell configuration.
#              Copy this file to your Pipeline_<MARKER>/ directory, rename it
#              config.sh, and fill in all PLACEHOLDER values.
#              Sourced automatically by main.sh and each pipeline script.
# Usage     : cp config/config_template.sh Pipeline_MYMARKER/config.sh
#             # then edit the copy — do NOT edit this template directly
# =============================================================================

# --- Project root (do not modify — resolved dynamically) ---------------------
# Always points to the MetaBarFlow/ root regardless of where the script is run
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# =============================================================================
# SECTION 1 — Marker identity
# =============================================================================

# Short marker label used to name output directories (e.g. 18S, 16S, ITS, COI)
MARKER="MARKER_PLACEHOLDER"

# Free-text description of the target group — used only in log headers
TAXONOMIC_GROUP="your_taxonomic_group"

# =============================================================================
# SECTION 2 — File naming
# =============================================================================

# Suffixes that identify forward and reverse read files.
# Change only if your sequencer uses a non-standard naming convention.
FASTQ_SUFFIX_R1="_R1.fastq.gz"
FASTQ_SUFFIX_R2="_R2.fastq.gz"

# =============================================================================
# SECTION 3 — Raw data locations
# =============================================================================

# Root folder containing all sequencing run directories for this marker.
# By convention: RAW_dataset/<MARKER>/
RAW_DIR="${PROJECT_ROOT}/RAW_dataset/${MARKER}"

# List every sequencing run directory that belongs to this marker.
# Each entry must contain paired FASTQ files (R1 + R2) directly inside it.
# Add or remove lines as needed; at least one entry is required.
RUN_DIRS=(
    "${RAW_DIR}/RunDir_1"   # REPLACE with actual run directory name
    # "${RAW_DIR}/RunDir_2" # uncomment and repeat for additional runs
)

# =============================================================================
# SECTION 4 — Pipeline output directories (do not modify)
# =============================================================================

PIPELINE_DIR="${PROJECT_ROOT}/results/${MARKER}"
SCRIPTS_DIR="${PROJECT_ROOT}/scripts"

STEP1_DIR="${PIPELINE_DIR}/01_merged"
STEP2_DIR="${PIPELINE_DIR}/02_QC"
STEP3_DIR="${PIPELINE_DIR}/03_trimmed"
STEP4_DIR="${PIPELINE_DIR}/04_ASV"
STEP5_DIR="${PIPELINE_DIR}/05_taxonomy"
LOG_DIR="${PIPELINE_DIR}/logs"
READ_TRACKING="${LOG_DIR}/read_tracking.tsv"

# =============================================================================
# SECTION 5 — Primers
# =============================================================================
# Provide exact primer sequences (5'→3', IUPAC ambiguity codes accepted).
# Reverse complements must be computed manually and checked — errors here
# will cause cutadapt to trim nothing or trim incorrectly.
#
# Pre-computed for common markers:
#   18S (E572F/E1009R) : FWD=CYGCGGTAATTCCAGCTC  REV=AYGGTATCTRATCRTCTTYG
#   ITS (ITS86F/ITS4R) : FWD=GTGAATCATCGAATCTTTGAA  REV=TCCTCCGCTTATTGATATGC
#   16S (515FB/926R)   : FWD=GTGYCAGCMGCCGCGGTAA  REV=CCGYCAATTYMTTTRAGTTT

PRIMER_FWD="FORWARD_PRIMER_SEQUENCE"
PRIMER_REV="REVERSE_PRIMER_SEQUENCE"
PRIMER_FWD_RC="FORWARD_PRIMER_REVERSE_COMPLEMENT"
PRIMER_REV_RC="REVERSE_PRIMER_REVERSE_COMPLEMENT"

# =============================================================================
# SECTION 6 — Cutadapt parameters
# =============================================================================

# Fraction of mismatches tolerated in primer matching (0.0–1.0).
# 0.2 = up to 20% errors — appropriate for degenerate/IUPAC primers.
CUTADAPT_ERROR_RATE=0.2

# Discard reads shorter than this after trimming (bp).
CUTADAPT_MIN_LENGTH=50

# Minimum overlap between primer and read for a match (bp).
# Lower values risk false positives; 10 is a safe default.
CUTADAPT_OVERLAP=10

# Number of CPU cores for cutadapt (set to 0 to use all available cores).
CUTADAPT_CORES=4

# Disallow insertions/deletions within primer region.
# true  = exact base-level matching within the primer (recommended for Illumina)
# false = allow indels (may be needed for long degenerate primers like 16S)
CUTADAPT_NO_INDELS=true

# How many times to search for adapter/primer per read.
# 1 = standard amplicon sequencing (single primer pair, no read-through)
# 2 = when read-through is expected (short amplicons on long reads, e.g. 16S V4-V5)
CUTADAPT_TIMES=1

# Set to true if the amplicon is shorter than 2× read length (read-through).
# When true, the pipeline passes both the forward and RC primer to cutadapt
# for linked adapter trimming. Typical cases:
#   18S V4  (~450 bp on 2×301): false
#   ITS2    (~411 bp on 2×301): true
#   16S V4-V5 (~267 bp on 2×301): true
CUTADAPT_READ_THROUGH=false

# =============================================================================
# SECTION 7 — FastQC / MultiQC
# =============================================================================

# Number of threads for FastQC (one thread per file up to this limit).
FASTQC_THREADS=4

# =============================================================================
# SECTION 8 — DADA2 parameters
# =============================================================================
# These values should be confirmed or adjusted after reviewing MultiQC output.
# The values below are conservative starting points.

# Truncation positions for forward and reverse reads (bp).
# Reads shorter than truncLen are discarded; longer reads are cut at this point.
# Set to 0 to DISABLE truncation (required for variable-length amplicons like ITS).
# Rule of thumb: keep positions where median quality ≥ 30 across all samples.
DADA2_TRUNC_LEN_FWD=230
DADA2_TRUNC_LEN_REV=180

# Maximum expected errors per read after truncation.
# Lower = stricter quality filter; 2 is standard for Illumina 2×300 data.
DADA2_MAX_EE_FWD=2
DADA2_MAX_EE_REV=2

# Truncate reads at the first base with quality score below this threshold.
# 2 is the DADA2 default (removes the lowest Phred scores only).
DADA2_TRUNC_Q=2

# Chimera removal method: "consensus" (default), "pooled", or "per-sample".
# "consensus" is fastest and works well for most datasets.
DADA2_CHIMERA_METHOD="consensus"

# Number of CPU threads for DADA2 (filterAndTrim, dada, mergePairs).
DADA2_THREADS=4

# =============================================================================
# SECTION 9 — Taxonomic assignment database
# =============================================================================
# Path to the reference database FASTA used by vsearch in step 5.
# Common databases and their download locations:
#   PR2 v5.1.0  (18S)  : https://github.com/pr2database/pr2database/releases
#   UNITE v9    (ITS)  : https://unite.ut.ee/repository.php
#   SILVA 138.2 (16S)  : https://www.arb-silva.de/download/arb-files/
# Place downloaded files in: databases/<DB_NAME>/

DB_PATH="${PROJECT_ROOT}/databases/YOUR_DB/your_database_file.fasta"

# Minimum sequence identity for a VSEARCH hit (0.0–1.0).
# 0.80 = 80% identity — appropriate for species-level assignment.
VSEARCH_IDENTITY=0.80

# Number of CPU threads for VSEARCH.
VSEARCH_THREADS=4

# Maximum number of accepted hits per query sequence.
# 1 returns only the best hit.
VSEARCH_MAXACCEPTS=1

# Maximum number of rejected hits before giving up on a query.
# 32 is the VSEARCH default and balances speed vs. sensitivity.
VSEARCH_MAXREJECTS=32

# =============================================================================
# SECTION 10 — Create pipeline directories (do not modify)
# =============================================================================
mkdir -p "${STEP1_DIR}" "${STEP2_DIR}" "${STEP3_DIR}" \
         "${STEP4_DIR}" "${STEP5_DIR}" "${LOG_DIR}"
