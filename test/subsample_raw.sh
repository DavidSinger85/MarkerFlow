#!/usr/bin/env bash
# =============================================================================
# Script    : subsample_raw.sh
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Description: Shrink the bundled example FASTQ files under RAW_dataset/ to a
#              fixed maximum number of read pairs each, so the repository stays
#              small enough for GitHub while still running the full pipeline
#              end to end. This documents exactly how the bundled example data
#              was produced from the original full-size sequencing files.
#
#              Real sequencing data (not synthetic). Each file is truncated to
#              its first CAP reads; because R1 and R2 are stored in the same
#              order, the first CAP reads of each are the same read pairs, so
#              pairing is preserved. Files already at or below CAP are left
#              untouched. Re-running is a no-op.
#
# Usage     : bash test/subsample_raw.sh [CAP] [ROOT]
#             CAP  = max read pairs to keep per file   (default 3000)
#             ROOT = directory to process in place     (default RAW_dataset)
# Example   : bash test/subsample_raw.sh 3000 RAW_dataset
# =============================================================================

set -eu

CAP="${1:-3000}"
ROOT="${2:-RAW_dataset}"

if [[ ! -d "${ROOT}" ]]; then
    echo "ERROR: directory not found: ${ROOT}" >&2
    echo "Run from the MetaBarFlow project root." >&2
    exit 1
fi

echo "============================================================================="
echo " Subsampling FASTQ under '${ROOT}' to a cap of ${CAP} read pairs per file"
echo "============================================================================="

lines_cap=$((CAP * 4))
total_before=0
total_after=0
n_trimmed=0
n_kept=0

# Iterate over every gzipped FASTQ, sorted for stable output.
while IFS= read -r -d '' f; do
    # Read count = lines / 4. Note: no `set -o pipefail` here on purpose —
    # `head` closing the pipe makes zcat exit 141, which is harmless.
    reads_before=$(( $(zcat "${f}" | wc -l) / 4 ))
    total_before=$((total_before + reads_before))

    if (( reads_before > CAP )); then
        tmp="${f}.tmp.$$"
        zcat "${f}" | head -n "${lines_cap}" | gzip > "${tmp}"
        mv "${tmp}" "${f}"
        reads_after=${CAP}
        n_trimmed=$((n_trimmed + 1))
        status="trimmed"
    else
        reads_after=${reads_before}
        n_kept=$((n_kept + 1))
        status="kept   "
    fi
    total_after=$((total_after + reads_after))

    printf "  [%s] %7d -> %7d reads   %s\n" \
        "${status}" "${reads_before}" "${reads_after}" "${f}"
done < <(find "${ROOT}" -name '*.fastq.gz' -print0 | sort -z)

echo "-----------------------------------------------------------------------------"
printf "  Files trimmed : %d\n" "${n_trimmed}"
printf "  Files kept    : %d\n" "${n_kept}"
printf "  Reads  : %d -> %d\n" "${total_before}" "${total_after}"
printf "  Size   : %s\n" "$(du -sh "${ROOT}" | cut -f1)"
echo "============================================================================="
