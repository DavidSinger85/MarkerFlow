#!/usr/bin/env bash
# =============================================================================
# Script    : 03_trim_primers.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Description: Trims primers from paired-end FastQ.gz files using Cutadapt.
#              Parameters fully driven by config.sh (marker-specific).
#              Handles read-through trimming when detected (e.g. ITS).
#              Discards read pairs where primers are not detected.
#              Runs post-trim QC at the end via 02_quality_check.sh.
# Usage     : bash scripts/03_trim_primers.sh <path/to/config.sh>
# Dependencies: cutadapt (>= 3.5)
# =============================================================================

set -euo pipefail

# --- Load configuration ------------------------------------------------------
CONFIG="${1:?ERROR: config.sh path required. Usage: bash 03_trim_primers.sh <path/to/config.sh>}"
CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"
source "${CONFIG}"

# --- Load shared utilities ---------------------------------------------------
source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# --- Log setup ---------------------------------------------------------------
PIPELINE_LOGFILE="${LOG_DIR}/03_trim_primers.log"
REPORTFILE="${LOG_DIR}/03_trim_primers_report.tsv"
log_header "03_trim_primers.sh" "${MARKER}" "${STEP1_DIR}" "${STEP3_DIR}"
log "Primer FWD      : ${PRIMER_FWD}"
log "Primer REV      : ${PRIMER_REV}"
log "Error rate      : ${CUTADAPT_ERROR_RATE}"
log "Min length      : ${CUTADAPT_MIN_LENGTH}"
log "Times           : ${CUTADAPT_TIMES}"
log "No-indels       : ${CUTADAPT_NO_INDELS}"
log "Read-through    : ${CUTADAPT_READ_THROUGH}"

# --- Dependency check --------------------------------------------------------
check_dependencies cutadapt

# --- Input check -------------------------------------------------------------
R1_FILES=("${STEP1_DIR}"/*"${FASTQ_SUFFIX_R1}")
if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    log "ERROR: no R1 FastQ.gz files found in ${STEP1_DIR} — aborting."
    exit 1
fi
log "Samples found   : ${#R1_FILES[@]}"

# --- Build Cutadapt argument string ------------------------------------------
build_cutadapt_args() {
    local ARGS=""
    ARGS+=" -g ${PRIMER_FWD}"
    ARGS+=" -G ${PRIMER_REV}"
    if [[ "${CUTADAPT_READ_THROUGH}" == "true" ]]; then
        ARGS+=" -a ${PRIMER_REV_RC}"
        ARGS+=" -A ${PRIMER_FWD_RC}"
    fi
    ARGS+=" --times ${CUTADAPT_TIMES}"
    ARGS+=" --error-rate ${CUTADAPT_ERROR_RATE}"
    ARGS+=" --overlap ${CUTADAPT_OVERLAP}"
    ARGS+=" --minimum-length ${CUTADAPT_MIN_LENGTH}"
    ARGS+=" --cores ${CUTADAPT_CORES}"
    ARGS+=" --discard-untrimmed"
    [[ "${CUTADAPT_NO_INDELS}" == "true" ]] && ARGS+=" --no-indels"
    echo "${ARGS}"
}

CUTADAPT_ARGS=$(build_cutadapt_args)
log "Cutadapt args   :${CUTADAPT_ARGS}"
log ""

# --- Report header -----------------------------------------------------------
echo -e "sample\treads_input\treads_trimmed\treads_discarded\tpct_retained" > "${REPORTFILE}"

# --- Trim primers ------------------------------------------------------------
log "----- Running Cutadapt : $(date) -----"

TOTAL_IN=0
TOTAL_OUT=0

for R1 in "${R1_FILES[@]}"; do

    R2="${R1/${FASTQ_SUFFIX_R1}/${FASTQ_SUFFIX_R2}}"
    if [[ ! -f "${R2}" ]]; then
        log "WARNING: R2 not found for $(basename "${R1}") — skipping."
        continue
    fi

    SAMPLE=$(basename "${R1}" "${FASTQ_SUFFIX_R1}")
    OUT_R1="${STEP3_DIR}/${SAMPLE}_R1.fastq.gz"
    OUT_R2="${STEP3_DIR}/${SAMPLE}_R2.fastq.gz"
    TMP_LOG="${LOG_DIR}/.cutadapt_${SAMPLE}.log"

    cutadapt \
        ${CUTADAPT_ARGS} \
        --output        "${OUT_R1}" \
        --paired-output "${OUT_R2}" \
        "${R1}" "${R2}" \
        > "${TMP_LOG}" 2>&1

    READS_IN=$(grep  "Total read pairs processed:" "${TMP_LOG}" | grep -oP '[\d,]+' | tr -d ',')
    READS_OUT=$(grep "Pairs written (passing filters):" "${TMP_LOG}" | grep -oP '[\d,]+' | head -1 | tr -d ',')
    READS_DISC=$(( READS_IN - READS_OUT ))
    PCT=$(awk "BEGIN {printf \"%.1f\", (${READS_OUT}/${READS_IN})*100}")

    echo -e "${SAMPLE}\t${READS_IN}\t${READS_OUT}\t${READS_DISC}\t${PCT}%" >> "${REPORTFILE}"
    log "  ${SAMPLE}: ${READS_IN} in → ${READS_OUT} retained (${PCT}%)"

    TOTAL_IN=$(( TOTAL_IN + READS_IN ))
    TOTAL_OUT=$(( TOTAL_OUT + READS_OUT ))

    update_read_tracking "${READ_TRACKING}" "${SAMPLE}" "step03_trimmed" "${READS_OUT}"

    cat "${TMP_LOG}" >> "${PIPELINE_LOGFILE}"
    rm -f "${TMP_LOG}"

done

# --- Global summary ----------------------------------------------------------
TOTAL_DISC=$(( TOTAL_IN - TOTAL_OUT ))
PCT_GLOBAL=$(awk "BEGIN {printf \"%.1f\", (${TOTAL_OUT}/${TOTAL_IN})*100}")

log ""
log "===== Trimming summary ==============================================="
log "  Total reads input     : ${TOTAL_IN}"
log "  Total reads retained  : ${TOTAL_OUT} (${PCT_GLOBAL}%)"
log "  Total reads discarded : ${TOTAL_DISC}"
log "======================================================================"
log_footer "03_trim_primers.sh" "${REPORTFILE}"

# --- Post-trim quality check -------------------------------------------------
log "----- Running post-trim QC : $(date) -----"
bash "$(dirname "${BASH_SOURCE[0]}")/02_quality_check.sh" "${CONFIG}" "${STEP3_DIR}"