#!/usr/bin/env bash
# =============================================================================
# Script    : 01_merge_fastq.sh
# Author    : David Singer
# Date      : 2026-03-01
# Description: Merges paired-end FastQ.gz files across N sequencing runs
#              defined in config.sh (RUN_DIRS array). Sample names are
#              auto-detected by stripping the standard Illumina suffix
#              ({sample}_S{n}_L{lane}_R{1/2}_{set}.fastq.gz) and the longest
#              common prefix across all files. Outputs standardized filenames:
#              {sample}_R1.fastq.gz / {sample}_R2.fastq.gz. Includes integrity
#              checks (gzip validation, read count verification, duplicate
#              detection).
# Usage     : bash scripts/01_merge_fastq.sh <path/to/config.sh>
# =============================================================================

set -euo pipefail

# --- Load configuration ------------------------------------------------------
CONFIG="${1:?ERROR: config.sh path required. Usage: bash 01_merge_fastq.sh <path/to/config.sh>}"
CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"
source "${CONFIG}"

# --- Load shared utilities ---------------------------------------------------
source "$(dirname "${BASH_SOURCE[0]}")/utils.sh"

# --- Log setup ---------------------------------------------------------------
PIPELINE_LOGFILE="${LOG_DIR}/01_merge_fastq.log"
REPORTFILE="${LOG_DIR}/01_merge_fastq_report.tsv"
log_header "01_merge_fastq.sh" "${MARKER}" "${RAW_DIR}" "${STEP1_DIR}"
log "Runs          : ${#RUN_DIRS[@]}"
for i in "${!RUN_DIRS[@]}"; do
    log "  Run $((i+1))      : ${RUN_DIRS[$i]}"
done

# --- Counters ----------------------------------------------------------------
COUNT_MERGED=0
COUNT_UNIQUE=0
COUNT_ERRORS=0

# --- Report header -----------------------------------------------------------
echo -e "output_file\truns_merged\tread_count_check\tstatus" > "${REPORTFILE}"

# =============================================================================
# Helper: strip sequencer suffix from a filename to get the bare sample token.
# Handles standard Illumina: {sample}_S{n}_L{lane}_R{1/2}_{set}.fastq.gz
# Falls back to stripping _R{1/2}[_set].fastq.gz for simpler conventions.
# =============================================================================
strip_suffix() {
    local BASE
    BASE=$(basename "$1")
    local STRIPPED
    STRIPPED=$(echo "${BASE}" | sed -E 's/_S[0-9]+_L[0-9]+_R[12][_0-9]*\.fastq\.gz$//')
    if [[ "${STRIPPED}" == "${BASE}" ]]; then
        STRIPPED=$(echo "${BASE}" | sed -E 's/_R[12][_0-9]*\.fastq\.gz$//')
    fi
    echo "${STRIPPED}"
}

# =============================================================================
# Helper: longest common prefix of all arguments
# =============================================================================
lcp() {
    local LCP="${1}"
    shift
    for STR in "$@"; do
        while [[ "${STR}" != "${LCP}"* ]]; do
            LCP="${LCP%?}"
            [[ -z "${LCP}" ]] && echo "" && return
        done
    done
    echo "${LCP}"
}

# =============================================================================
# Helper: extract read direction (R1 or R2) from filename
# =============================================================================
extract_read() {
    echo "$1" | grep -oP 'R[12]'
}

# =============================================================================
# Auto-detect common prefix across all runs
# Strips the Illumina suffix from every filename, computes the longest common
# prefix, then trims it back to the last word separator (_ - .) so the
# remaining token is the clean sample identifier.
# =============================================================================
log "----- Auto-detecting sample names -----"

ALL_STRIPPED=()
for _RUN_DIR in "${RUN_DIRS[@]}"; do
    [[ -d "${_RUN_DIR}" ]] || continue
    for _FILE in "${_RUN_DIR}"/*.fastq.gz; do
        [[ -f "${_FILE}" ]] || continue
        ALL_STRIPPED+=("$(strip_suffix "${_FILE}")")
    done
done

if [[ ${#ALL_STRIPPED[@]} -eq 0 ]]; then
    log "ERROR: No FastQ.gz files found in any run directory — aborting."
    exit 1
fi

# Deduplicate (R1 and R2 of the same sample produce the same stripped name)
mapfile -t UNIQUE_STRIPPED < <(printf '%s\n' "${ALL_STRIPPED[@]}" | sort -u)

# Compute LCP then trim to last word separator so we don't cut mid-token
COMMON_PREFIX=$(lcp "${UNIQUE_STRIPPED[@]}")
if [[ -n "${COMMON_PREFIX}" ]]; then
    COMMON_PREFIX=$(echo "${COMMON_PREFIX}" | sed -E 's/[^_.\-]+$//')
fi
log "  Common prefix stripped : '${COMMON_PREFIX}'"
log "  Sample IDs detected    : ${#UNIQUE_STRIPPED[@]}"
for S in "${UNIQUE_STRIPPED[@]}"; do
    log "    ${S#${COMMON_PREFIX}}"
done

# =============================================================================
# Step 1: Index all runs — build a 2D index: IDX[RUN_NUM:SAMPLE_READ] → path
#         Detect duplicates within the same run
# =============================================================================
declare -A IDX
declare -A ALL_KEYS

for RUN_I in "${!RUN_DIRS[@]}"; do
    RUN_DIR="${RUN_DIRS[$RUN_I]}"
    if [[ ! -d "${RUN_DIR}" ]]; then
        log "ERROR: run directory not found: ${RUN_DIR} — aborting."
        exit 1
    fi
    for FILE in "${RUN_DIR}"/*.fastq.gz; do
        BASE=$(basename "${FILE}")
        SAMPLE=$(strip_suffix "${FILE}")
        SAMPLE="${SAMPLE#${COMMON_PREFIX}}"
        READ=$(extract_read "${BASE}")
        if [[ -z "${SAMPLE}" || -z "${READ}" ]]; then
            log "WARNING: cannot parse filename in run $((RUN_I+1)): ${BASE} — skipping."
            continue
        fi
        KEY="${SAMPLE}_${READ}"
        IDXKEY="${RUN_I}:${KEY}"
        if [[ -v IDX["${IDXKEY}"] ]]; then
            log "ERROR: duplicate sample+read in run $((RUN_I+1)): ${BASE} — aborting."
            exit 1
        fi
        IDX["${IDXKEY}"]="${FILE}"
        ALL_KEYS["${KEY}"]=1
    done
done

# =============================================================================
# Step 2: For each unique sample+read, merge across all runs that have it
# =============================================================================
TMPDIR_MERGE="${STEP1_DIR}/.tmp_merge"
mkdir -p "${TMPDIR_MERGE}"

for KEY in "${!ALL_KEYS[@]}"; do
    SAMPLE="${KEY%_*}"
    READ="${KEY##*_}"
    OUTFILE="${STEP1_DIR}/${SAMPLE}_${READ}.fastq.gz"

    FILES_TO_MERGE=()
    RUNS_PRESENT=()
    for RUN_I in "${!RUN_DIRS[@]}"; do
        IDXKEY="${RUN_I}:${KEY}"
        if [[ -v IDX["${IDXKEY}"] ]]; then
            FILES_TO_MERGE+=("${IDX[${IDXKEY}]}")
            RUNS_PRESENT+=("$((RUN_I+1))")
        fi
    done

    RUNS_LABEL=$(IFS=+; echo "${RUNS_PRESENT[*]}")

    if [[ ${#FILES_TO_MERGE[@]} -eq 1 ]]; then
        cp "${FILES_TO_MERGE[0]}" "${OUTFILE}"
        check_gzip "${OUTFILE}"
        RC_FINAL=$(count_reads "${OUTFILE}")
        echo -e "${SAMPLE}_${READ}.fastq.gz\trun ${RUNS_LABEL} only\tOK (${RC_FINAL} reads)\tOK" >> "${REPORTFILE}"
        COUNT_UNIQUE=$(( COUNT_UNIQUE + 1 ))
        [[ "${READ}" == "R1" ]] && update_read_tracking "${READ_TRACKING}" "${SAMPLE}" "step01_merged" "${RC_FINAL}"
    else
        EXPECTED=0
        for F in "${FILES_TO_MERGE[@]}"; do
            RC=$(count_reads "${F}")
            EXPECTED=$(( EXPECTED + RC ))
        done

        cat "${FILES_TO_MERGE[@]}" > "${OUTFILE}"
        check_gzip "${OUTFILE}"
        RC_FINAL=$(count_reads "${OUTFILE}")

        if [[ "${RC_FINAL}" -ne "${EXPECTED}" ]]; then
            log "ERROR: read count mismatch for ${SAMPLE}_${READ}: expected ${EXPECTED}, got ${RC_FINAL}"
            echo -e "${SAMPLE}_${READ}.fastq.gz\t${RUNS_LABEL}\tFAIL (expected ${EXPECTED}, got ${RC_FINAL})\tERROR" >> "${REPORTFILE}"
            COUNT_ERRORS=$(( COUNT_ERRORS + 1 ))
        else
            echo -e "${SAMPLE}_${READ}.fastq.gz\t${RUNS_LABEL}\tOK (${RC_FINAL} reads)\tOK" >> "${REPORTFILE}"
            COUNT_MERGED=$(( COUNT_MERGED + 1 ))
            [[ "${READ}" == "R1" ]] && update_read_tracking "${READ_TRACKING}" "${SAMPLE}" "step01_merged" "${RC_FINAL}"
        fi
    fi
done

# =============================================================================
# Step 3: Cleanup
# =============================================================================
rm -rf "${TMPDIR_MERGE}"

# =============================================================================
# Step 4: Summary
# =============================================================================
log ""
log "===== Merge summary =================================================="
log "  Samples merged (2+ runs)   : ${COUNT_MERGED}"
log "  Samples unique (1 run)     : ${COUNT_UNIQUE}"
log "  Errors detected            : ${COUNT_ERRORS}"
log "======================================================================"
log_footer "01_merge_fastq.sh" "${REPORTFILE}"
