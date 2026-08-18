#!/usr/bin/env bash
# =============================================================================
# File      : utils.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Description: Shared utility functions for the metabarcoding pipeline.
#              Uses log() wrapper instead of exec redirection for WSL2
#              compatibility and safe nested script calls.
# Usage     : source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"
# =============================================================================

# =============================================================================
# PIPELINE_LOGFILE
# Set by each script before calling any log function:
#   PIPELINE_LOGFILE="${LOG_DIR}/02_quality_check.log"
# =============================================================================
PIPELINE_LOGFILE=""

# =============================================================================
# log
# Writes to terminal AND PIPELINE_LOGFILE. No exec redirection — WSL2 safe.
# =============================================================================
log() {
    echo "$@"
    [[ -n "${PIPELINE_LOGFILE}" ]] && echo "$@" >> "${PIPELINE_LOGFILE}"
}

# =============================================================================
# log_header
# Arguments: $1 script name  $2 marker  $3 input dir  $4 output dir
# =============================================================================
log_header() {
    local SCRIPT="$1" MARKER="$2" INPUT="$3" OUTPUT="$4"
    log "===== ${SCRIPT} started : $(date) ====="
    log "Marker        : ${MARKER}"
    log "Input folder  : ${INPUT}"
    log "Output folder : ${OUTPUT}"
}

# =============================================================================
# log_footer
# Arguments: $1 script name  $2 report file (optional)
# =============================================================================
log_footer() {
    local SCRIPT="$1" REPORTFILE="${2:-}"
    [[ -n "${REPORTFILE}" ]] && log "  Report : ${REPORTFILE}"
    log "  Log    : ${PIPELINE_LOGFILE}"
    log "===== ${SCRIPT} finished : $(date) ====="
}

# =============================================================================
# check_dependencies
# Arguments: $@ list of tool names
# =============================================================================
check_dependencies() {
    local MISSING=0
    log "--- Dependency check ---"
    for CMD in "$@"; do
        if ! command -v "${CMD}" &>/dev/null; then
            log "ERROR: ${CMD} not found in PATH — aborting."
            MISSING=$(( MISSING + 1 ))
        else
            # || true: prevents set -e on non-zero exit (e.g. fastqc v0.11.x)
            local VERSION
            VERSION=$(${CMD} --version 2>&1 | head -1 || true)
            log "  ${CMD} : ${VERSION}"
        fi
    done
    if [[ "${MISSING}" -gt 0 ]]; then exit 1; fi
}

# =============================================================================
# check_gzip
# Arguments: $1 path to .gz file
# =============================================================================
check_gzip() {
    if ! gzip -t "$1" 2>/dev/null; then
        log "ERROR: corrupted gzip file: $1"
        return 1
    fi
}

# =============================================================================
# count_reads
# Arguments: $1 path to .fastq.gz file
# =============================================================================
count_reads() {
    zcat "$1" | awk 'END {print NR/4}'
}

# =============================================================================
# update_read_tracking
# Arguments: $1 tracking file  $2 sample  $3 step name  $4 count
# =============================================================================
update_read_tracking() {
    local TRACKING_FILE="$1" SAMPLE="$2" STEP="$3" COUNT="$4"
    local TMP="${TRACKING_FILE}.tmp"

    if [[ ! -f "${TRACKING_FILE}" ]]; then
        echo -e "sample\t${STEP}" > "${TRACKING_FILE}"
    fi

    if ! head -1 "${TRACKING_FILE}" | grep -qP "\t${STEP}(\t|$)"; then
        awk -v step="${STEP}" '
            NR==1 { print $0 "\t" step; next }
                  { print $0 "\tNA" }
        ' "${TRACKING_FILE}" > "${TMP}" && mv "${TMP}" "${TRACKING_FILE}"
    fi

    if grep -qP "^${SAMPLE}\t" "${TRACKING_FILE}"; then
        awk -v sample="${SAMPLE}" -v step="${STEP}" -v count="${COUNT}" '
            NR==1 { for (i=1; i<=NF; i++) if ($i==step) col=i; print; next }
            $1==sample { $col=count; print; next }
            { print }
        ' OFS='\t' "${TRACKING_FILE}" > "${TMP}" && mv "${TMP}" "${TRACKING_FILE}"
    else
        local NCOLS COL_IDX NEW_ROW
        NCOLS=$(head -1 "${TRACKING_FILE}" | awk '{print NF}')
        COL_IDX=$(head -1 "${TRACKING_FILE}" | tr '\t' '\n' | grep -n "^${STEP}$" | cut -d: -f1)
        NEW_ROW="${SAMPLE}"
        for (( i=2; i<=NCOLS; i++ )); do
            [[ "${i}" -eq "${COL_IDX}" ]] && NEW_ROW+="\t${COUNT}" || NEW_ROW+="\tNA"
        done
        echo -e "${NEW_ROW}" >> "${TRACKING_FILE}"
    fi
}
