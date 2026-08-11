#!/usr/bin/env bash
# =============================================================================
# Script    : main.sh
# Author    : David Singer
# Date      : 2026-03-01
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Description: Master launcher for the metabarcoding pipeline. Executes one
#              or more steps sequentially for a given marker configuration.
#              Steps can be selected individually or combined to allow
#              parameter testing and re-runs of specific steps.
# Usage     : bash main.sh <path/to/config.sh> --steps <1,2,3,4,5>
# Examples  : bash main.sh Pipeline_18S/config.sh --steps 1,2,3,4,5
#             bash main.sh Pipeline_18S/config.sh --steps 3
#             bash main.sh Pipeline_ITS/config.sh --steps 2,3
# =============================================================================

set -euo pipefail

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: bash main.sh <path/to/config.sh> --steps <1,2,3,4,5>"
    echo "  Steps:"
    echo "    1 — Merge FastQ files"
    echo "    2 — Quality check (FastQC / MultiQC)"
    echo "    3 — Primer trimming (Cutadapt)"
    echo "    4 — ASV generation (DADA2 via Rscript)"
    echo "    5 — Taxonomic assignment (VSEARCH)"
    echo "  Examples:"
    echo "    bash main.sh Pipeline_18S/config.sh --steps 1,2,3,4,5"
    echo "    bash main.sh Pipeline_18S/config.sh --steps 3"
    echo "    bash main.sh Pipeline_18S/config.sh --steps 2,3"
    exit 1
}

# --- Argument parsing --------------------------------------------------------
CONFIG="${1:-}"
STEPS_ARG="${3:-}"

[[ -z "${CONFIG}" || "${2:-}" != "--steps" || -z "${STEPS_ARG}" ]] && usage
[[ ! -f "${CONFIG}" ]] && { echo "ERROR: config not found: ${CONFIG}"; exit 1; }

CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"
source "${CONFIG}"

# --- Resolve scripts directory -----------------------------------------------
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Load shared utilities ---------------------------------------------------
source "${SCRIPTS_DIR}/utils.sh"

# --- Log setup ---------------------------------------------------------------
PIPELINE_LOGFILE="${LOG_DIR}/main.log"
log "===== main.sh started : $(date) ====="
log "Marker  : ${MARKER}"
log "Config  : ${CONFIG}"
log "Steps   : ${STEPS_ARG}"

# --- Parse steps into array --------------------------------------------------
IFS=',' read -ra STEPS <<< "${STEPS_ARG}"

for STEP in "${STEPS[@]}"; do
    if [[ ! "${STEP}" =~ ^[1-5]$ ]]; then
        log "ERROR: invalid step '${STEP}' — must be 1, 2, 3, 4, or 5."
        usage
    fi
done

# --- Helper: check if a step is requested ------------------------------------
run_step() {
    local S="$1"
    for STEP in "${STEPS[@]}"; do
        [[ "${STEP}" == "${S}" ]] && return 0
    done
    return 1
}

# --- Step 1: Merge FastQ -----------------------------------------------------
if run_step 1; then
    log "----- Step 1: Merge FastQ : $(date) -----"
    bash "${SCRIPTS_DIR}/01_merge_fastq.sh" "${CONFIG}"
    log "----- Step 1 complete -----"
fi

# --- Step 2: Quality check ---------------------------------------------------
if run_step 2; then
    log "----- Step 2: Quality check : $(date) -----"
    bash "${SCRIPTS_DIR}/02_quality_check.sh" "${CONFIG}"
    log "----- Step 2 complete -----"
fi

# --- Step 3: Primer trimming -------------------------------------------------
if run_step 3; then
    log "----- Step 3: Primer trimming : $(date) -----"
    bash "${SCRIPTS_DIR}/03_trim_primers.sh" "${CONFIG}"
    log "----- Step 3 complete -----"
fi

# --- Step 4: ASV generation (DADA2) ------------------------------------------
if run_step 4; then
    log "----- Step 4: ASV generation (DADA2) : $(date) -----"
    if ! command -v Rscript &>/dev/null; then
        log "ERROR: Rscript not found in PATH — aborting step 4."
        exit 1
    fi
    CONFIG_R="$(dirname "${CONFIG}")/config.R"
    [[ ! -f "${CONFIG_R}" ]] && { log "ERROR: config.R not found: ${CONFIG_R}"; exit 1; }
    Rscript "${SCRIPTS_DIR}/04_DADA2_ASV.R" "${CONFIG_R}"
    log "----- Step 4 complete -----"
fi

# --- Step 5: Taxonomic assignment --------------------------------------------
if run_step 5; then
    log "----- Step 5: Taxonomic assignment : $(date) -----"
    if ! command -v Rscript &>/dev/null; then
        log "ERROR: Rscript not found in PATH — aborting step 5."
        exit 1
    fi
    CONFIG_R="$(dirname "${CONFIG}")/config.R"
    [[ ! -f "${CONFIG_R}" ]] && { log "ERROR: config.R not found: ${CONFIG_R}"; exit 1; }
    Rscript "${SCRIPTS_DIR}/05_taxonomy.R" "${CONFIG_R}"
    log "----- Step 5 complete -----"
fi

# --- Done --------------------------------------------------------------------
log_footer "main.sh"