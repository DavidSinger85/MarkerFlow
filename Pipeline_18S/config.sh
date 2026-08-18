#!/usr/bin/env bash
# =============================================================================
# File      : config.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Description: Marker-specific configuration file for the 18S metabarcoding
#              pipeline (protists / testate amoebae). Sourced by all pipeline
#              scripts. To adapt for ITS or 16S, copy this file to the
#              corresponding Pipeline_ITS/ or Pipeline_16S/ directory and
#              update the parameters below.
# Usage     : source config.sh  (called automatically by each pipeline script)
# =============================================================================

# --- Project root (resolved dynamically — GitHub portable) -------------------
# Resolves to the parent directory of the Pipeline_18S/ folder
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Marker identity ---------------------------------------------------------
MARKER="18S"
TAXONOMIC_GROUP="protists_testate_amoebae"

# --- File naming -------------------------------------------------------------
FASTQ_SUFFIX_R1="_R1.fastq.gz"
FASTQ_SUFFIX_R2="_R2.fastq.gz"

# --- Raw data ----------------------------------------------------------------
RAW_DIR="${PROJECT_ROOT}/RAW_dataset/18S"

# Array of sequencing run directories — add or remove runs as needed
RUN_DIRS=(
    "${RAW_DIR}/Heger18SV4_1"
    "${RAW_DIR}/Heger18SV4_2"
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

# --- Primers (E572F / E1009R — Elwood et al. 1985, modified) ----------------
PRIMER_FWD="CYGCGGTAATTCCAGCTC"
PRIMER_REV="AYGGTATCTRATCRTCTTYG"
PRIMER_FWD_RC="GAGCTGGAATTACCGCRG"
PRIMER_REV_RC="CRAAGAYGATYAGATACCRT"

# --- Cutadapt parameters -----------------------------------------------------
CUTADAPT_ERROR_RATE=0.2
CUTADAPT_MIN_LENGTH=50
CUTADAPT_OVERLAP=10
CUTADAPT_CORES=4
CUTADAPT_NO_INDELS=true
CUTADAPT_TIMES=1
CUTADAPT_READ_THROUGH=false   # no read-through detected (amplicon ~450bp > 2x301bp)

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
DB_PATH="${PROJECT_ROOT}/databases/PR2/pr2_version_5.1.0_SSU_mothur.fasta"
VSEARCH_IDENTITY=0.80
VSEARCH_THREADS=4
VSEARCH_MAXACCEPTS=1
VSEARCH_MAXREJECTS=32

# --- Create all pipeline directories if not yet present ----------------------
mkdir -p "${STEP1_DIR}" "${STEP2_DIR}" "${STEP3_DIR}" \
         "${STEP4_DIR}" "${STEP5_DIR}" "${LOG_DIR}"