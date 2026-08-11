#!/usr/bin/env bash
# =============================================================================
# File      : config.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Description: Marker-specific configuration file for the 16S metabarcoding
#              pipeline (bacteria). Sourced by all pipeline scripts.
#              Targets the V4-V5 region using primers 515FB / 926R.
#              PROJECT_ROOT is resolved dynamically — fully portable.
# Usage     : source config.sh  (called automatically by each pipeline script)
# =============================================================================

# --- Project root (resolved dynamically — GitHub portable) -------------------
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Marker identity ---------------------------------------------------------
MARKER="16S"
TAXONOMIC_GROUP="bacteria"

# --- File naming -------------------------------------------------------------
FASTQ_SUFFIX_R1="_R1.fastq.gz"
FASTQ_SUFFIX_R2="_R2.fastq.gz"

# --- Raw data ----------------------------------------------------------------
RAW_DIR="${PROJECT_ROOT}/RAW_dataset/16S"

# Array of sequencing run directories — add or remove runs as needed
RUN_DIRS=(
    "${RAW_DIR}/HegerV4V5_1"
    "${RAW_DIR}/HegerV4V5_2"
)

# --- Pipeline directories ----------------------------------------------------
PIPELINE_DIR="${PROJECT_ROOT}/results/${MARKER}"
SCRIPTS_DIR="${PROJECT_ROOT}/scripts"

STEP1_DIR="${PIPELINE_DIR}/01_merged"
STEP2_DIR="${PIPELINE_DIR}/02_QC"
STEP3_DIR="${PIPELINE_DIR}/03_trimmed"
STEP4_DIR="${PIPELINE_DIR}/04_ASV"
STEP5_DIR="${PIPELINE_DIR}/05_taxonomy"
LOG_DIR="${PIPELINE_DIR}/logs"
READ_TRACKING="${LOG_DIR}/read_tracking.tsv"

# --- Primers (515FB / 926R — Parada et al. 2016 / Quince et al. 2011) -------
# Forward primer : 515FB
PRIMER_FWD="GTGYCAGCMGCCGCGGTAA"
# Reverse primer : 926R
PRIMER_REV="CCGYCAATTYMTTTRAGTTT"
# Reverse complements (verified manually)
PRIMER_FWD_RC="TTACCGCGGCKGCTGRCAC"
PRIMER_REV_RC="AAACTAAAKRAATTYRGGCGG"

# --- Cutadapt parameters -----------------------------------------------------
CUTADAPT_ERROR_RATE=0.2
CUTADAPT_MIN_LENGTH=50
CUTADAPT_OVERLAP=10
CUTADAPT_CORES=4
CUTADAPT_NO_INDELS=false
CUTADAPT_TIMES=2
CUTADAPT_READ_THROUGH=true    # read-through confirmed (~80%, amplicon ~267bp < 2x301bp)

# --- FastQC / MultiQC --------------------------------------------------------
FASTQC_THREADS=4

# --- DADA2 parameters (update after QC inspection) ---------------------------
DADA2_TRUNC_LEN_FWD=230
DADA2_TRUNC_LEN_REV=180
DADA2_MAX_EE_FWD=2
DADA2_MAX_EE_REV=2
DADA2_TRUNC_Q=2
DADA2_CHIMERA_METHOD="consensus"
DADA2_THREADS=4

# --- Taxonomic assignment (VSEARCH) ------------------------------------------
# Database path — adjust to local installation, do not push databases to GitHub
DB_PATH="${PROJECT_ROOT}/databases/SILVA/silva_nr99_v138.2_train_set.fa.gz"
VSEARCH_IDENTITY=0.80
VSEARCH_THREADS=4
VSEARCH_MAXACCEPTS=1
VSEARCH_MAXREJECTS=32

# --- Create all pipeline directories if not yet present ----------------------
mkdir -p "${STEP1_DIR}" "${STEP2_DIR}" "${STEP3_DIR}" \
         "${STEP4_DIR}" "${STEP5_DIR}" "${LOG_DIR}"