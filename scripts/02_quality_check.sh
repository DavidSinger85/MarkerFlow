#!/usr/bin/env bash
# =============================================================================
# Script    : 02_quality_check.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Description: Runs FastQC and MultiQC. Pre-trim (default) or post-trim mode.
#              Pre-trim output  : STEP2_DIR/pre_trim/
#              Post-trim output : STEP2_DIR/post_trim/
# Usage     : bash scripts/02_quality_check.sh <path/to/config.sh> [input_dir]
# Dependencies: fastqc, multiqc
# =============================================================================

set -euo pipefail

# --- Load configuration ------------------------------------------------------
CONFIG="${1:?ERROR: config.sh path required.}"
CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"
source "${CONFIG}"

# --- Load shared utilities ---------------------------------------------------
source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# --- Resolve input directory and QC mode -------------------------------------
if [[ -n "${2:-}" ]]; then
    INPUT_DIR="${2}"
    QC_MODE="post_trim"
else
    INPUT_DIR="${STEP1_DIR}"
    QC_MODE="pre_trim"
fi
OUTPUT_DIR="${STEP2_DIR}/${QC_MODE}"
mkdir -p "${OUTPUT_DIR}"

# --- Log setup ---------------------------------------------------------------
PIPELINE_LOGFILE="${LOG_DIR}/02_quality_check_${QC_MODE}.log"
log_header "02_quality_check.sh [${QC_MODE}]" "${MARKER}" "${INPUT_DIR}" "${OUTPUT_DIR}"

# --- Dependency check --------------------------------------------------------
check_dependencies fastqc multiqc

# --- Input check -------------------------------------------------------------
FASTQ_FILES=("${INPUT_DIR}"/*.fastq.gz)
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    log "ERROR: no FastQ.gz files found in ${INPUT_DIR} — aborting."
    exit 1
fi
log "Samples found : ${#FASTQ_FILES[@]} files"

# --- Run FastQC --------------------------------------------------------------
log "----- Running FastQC : $(date) -----"
fastqc \
    --outdir  "${OUTPUT_DIR}" \
    --threads "${FASTQC_THREADS}" \
    --quiet \
    "${FASTQ_FILES[@]}" || true
log "FastQC complete."

# --- Run MultiQC -------------------------------------------------------------
log "----- Running MultiQC : $(date) -----"
multiqc \
    "${OUTPUT_DIR}" \
    --outdir   "${OUTPUT_DIR}/multiqc" \
    --filename "multiqc_${MARKER}_${QC_MODE}_report" \
    --quiet \
    --force
log "MultiQC complete."
log "Report : ${OUTPUT_DIR}/multiqc/multiqc_${MARKER}_${QC_MODE}_report.html"

# --- Completion --------------------------------------------------------------
log_footer "02_quality_check.sh [${QC_MODE}]"