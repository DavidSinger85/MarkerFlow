#!/usr/bin/env bash
# =============================================================================
# File      : config.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Description: Marker-specific configuration file for the ITS metabarcoding
#              pipeline (fungi). Sourced by all pipeline scripts. Based on
#              primers ITS86F / ITS4R targeting the ITS2 region.
# Usage     : source config.sh  (called automatically by each pipeline script)
# =============================================================================

# --- Project root (resolved dynamically — GitHub portable) -------------------
# Resolves to the parent directory of the Pipeline_18S/ folder
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Marker identity ---------------------------------------------------------
MARKER="ITS"
TAXONOMIC_GROUP="fungi"

# --- File naming -------------------------------------------------------------
FASTQ_SUFFIX_R1="_R1.fastq.gz"
FASTQ_SUFFIX_R2="_R2.fastq.gz"

# --- Raw data ----------------------------------------------------------------
RAW_DIR="${PROJECT_ROOT}/RAW_dataset/ITS"

# Array of sequencing run directories — add or remove runs as needed
RUN_DIRS=(
    "${RAW_DIR}/Heger_ITS_1"
    "${RAW_DIR}/Heger_ITS_2"
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

# --- Primers (ITS86F / ITS4R — Turenne et al. 1999 / White et al. 1990) -----
# Forward primer : ITS86F
PRIMER_FWD="GTGAATCATCGAATCTTTGAA"
# Reverse primer : ITS4R
PRIMER_REV="TCCTCCGCTTATTGATATGC"
# Reverse complements (verified manually)
PRIMER_FWD_RC="TTCAAAGATTCGATGATTCAC"
PRIMER_REV_RC="GCATATCAATAAGCGGAGGA"

# --- Cutadapt parameters -----------------------------------------------------
CUTADAPT_ERROR_RATE=0.2
CUTADAPT_MIN_LENGTH=50
CUTADAPT_OVERLAP=10
CUTADAPT_CORES=4
CUTADAPT_NO_INDELS=true
CUTADAPT_TIMES=1
CUTADAPT_READ_THROUGH=true   # no read-through detected (amplicon ~411bp > 2x301bp)

# --- FastQC / MultiQC --------------------------------------------------------
FASTQC_THREADS=4

# --- DADA2 parameters (update after QC inspection) ---------------------------
# Note: ITS amplicon length is highly variable — truncation is typically
# disabled (set to 0) for ITS; use DADA2_MAX_EE filtering only
DADA2_TRUNC_LEN_FWD=0
DADA2_TRUNC_LEN_REV=0
DADA2_MAX_EE_FWD=2
DADA2_MAX_EE_REV=2
DADA2_TRUNC_Q=2
DADA2_CHIMERA_METHOD="consensus"
DADA2_THREADS=4

# --- Taxonomic assignment (VSEARCH) ------------------------------------------
DB_PATH="${PROJECT_ROOT}/databases/UNITE/sh_general_release_dynamic_19.02.2025.fasta"
VSEARCH_IDENTITY=0.80
VSEARCH_THREADS=4
VSEARCH_MAXACCEPTS=1
VSEARCH_MAXREJECTS=32

# --- Create all pipeline directories if not yet present ----------------------
mkdir -p "${STEP1_DIR}" "${STEP2_DIR}" "${STEP3_DIR}" \
         "${STEP4_DIR}" "${STEP5_DIR}" "${LOG_DIR}"