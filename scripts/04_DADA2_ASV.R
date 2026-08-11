# =============================================================================
# File      : 04_DADA2_ASV.R
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Author    : David Singer
# Date      : 2026-03-03
# Description: DADA2 ASV generation pipeline.
#              Reads trimmed FastQ from STEP3_DIR, runs filterAndTrim,
#              learnErrors, dada, mergePairs, makeSequenceTable,
#              removeBimeraDenovo. Outputs:
#                - [MARKER]_MR.csv       : community matrix (samples x ASVs)
#                - [MARKER]_Fasta.fasta  : ASV sequences
#                - seqtab_nochim.rds     : sorted/renamed seqtab for step 05
#              Updates read_tracking.tsv.
# Usage     : Rscript 04_DADA2_ASV.R pipeline_18S/config.R
# =============================================================================

# --- Dependencies ------------------------------------------------------------
suppressPackageStartupMessages({
    library(dada2)
})

# =============================================================================
# UTILITIES
# =============================================================================

#' Write sequences to FASTA file — base R, no seqinr dependency
#' @param sequences character vector of sequences
#' @param names     character vector of sequence names
#' @param file      output file path
write_fasta <- function(sequences, names, file) {
    con <- file(file, open = "w")
    for (i in seq_along(sequences)) {
        writeLines(paste0(">", names[i]), con)
        writeLines(sequences[i], con)
    }
    close(con)
}

# =============================================================================
# LOGGING
# =============================================================================

#' Log a message to terminal and log file via sink()
#' @param ... message components passed to paste0()
log_r <- function(...) {
    cat(paste0(..., "\n"))
}

# =============================================================================
# READ TRACKING
# =============================================================================

#' Update read_tracking.tsv with new columns
#' @param tracking_file path to read_tracking.tsv
#' @param new_data named data.frame with sample as rowname
update_read_tracking <- function(tracking_file, new_data) {
    if (file.exists(tracking_file)) {
        # Read first column as character to preserve leading zeros (e.g. "033")
        # row.names=1 silently converts "033" -> 33 -> "33", causing duplicates
        raw <- read.table(tracking_file,
                          header      = TRUE,
                          sep         = "\t",
                          check.names = FALSE,
                          colClasses  = "character",
                          row.names   = NULL)
        existing <- raw[, -1, drop = FALSE]
        rownames(existing) <- raw[, 1]
        for (col in colnames(new_data)) {
            existing[rownames(new_data), col] <- new_data[, col]
        }
        write.table(existing, tracking_file,
                    sep       = "\t",
                    quote     = FALSE,
                    col.names = NA)
    } else {
        write.table(new_data, tracking_file,
                    sep       = "\t",
                    quote     = FALSE,
                    col.names = NA)
    }
}

# =============================================================================
# INITIALISATION
# =============================================================================

# --- Parse config path argument ----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 04_DADA2_ASV.R <path/to/config.R>")
}
config_path  <- normalizePath(args[1], mustWork = TRUE)
source(config_path)

# Coerce POOLING to correct type — config.R may read logical as character
if (is.character(POOLING)) {
    POOLING <- switch(POOLING,
        "TRUE"   = TRUE,
        "FALSE"  = FALSE,
        "pseudo" = "pseudo",
        stop("Invalid POOLING value in config.R: must be TRUE, FALSE, or 'pseudo'")
    )
}

# Fallbacks for parameters added in later config versions
if (!exists("THREADS_LEARN"))       THREADS_LEARN       <- 1
if (!exists("MIN_OVERLAP"))         MIN_OVERLAP         <- 12
if (!exists("FASTQ_STRIP_PATTERN_R")) FASTQ_STRIP_PATTERN_R <- "_R1\\.fastq\\.gz$"

# --- Setup log file ----------------------------------------------------------
log_file <- file.path(LOG_DIR, "04_DADA2_ASV.log")
dir.create(LOG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(STEP4_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STEP5_DIR, recursive = TRUE, showWarnings = FALSE)

sink(file = log_file, append = FALSE, split = TRUE)  # split=TRUE → terminal + file simultaneously

log_r("=============================================================================")
log_r("  04_DADA2_ASV.R started : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_r("=============================================================================")
log_r("  Marker       : ", MARKER)
log_r("  Config       : ", config_path)
log_r("  Input        : ", STEP3_DIR)
log_r("  Output       : ", STEP4_DIR)
log_r("  Threads      : ", THREADS)
log_r("  Pooling      : ", POOLING)
log_r("  Chimera      : ", CHIMERA_METHOD)
log_r("-----------------------------------------------------------------------------")

# =============================================================================
# STEP 1 — COLLECT INPUT FILES
# =============================================================================

log_r("----- Step 1: Collect input files : ", format(Sys.time(), "%H:%M:%S"), " -----")

fastq_pattern_r2 <- sub("R1", "R2", FASTQ_STRIP_PATTERN_R, fixed = TRUE)
fwd_files <- sort(list.files(STEP3_DIR, pattern = FASTQ_STRIP_PATTERN_R, full.names = TRUE))
rev_files <- sort(list.files(STEP3_DIR, pattern = fastq_pattern_r2,      full.names = TRUE))

if (length(fwd_files) == 0) stop("No R1 FastQ files found in ", STEP3_DIR)
if (length(fwd_files) != length(rev_files)) stop("Mismatch between R1 and R2 file counts")

sample_names <- gsub(FASTQ_STRIP_PATTERN_R, "", basename(fwd_files))
log_r("  Samples found : ", length(sample_names))

# =============================================================================
# STEP 2 — FILTER AND TRIM
# =============================================================================

log_r("----- Step 2: filterAndTrim : ", format(Sys.time(), "%H:%M:%S"), " -----")
log_r("  TRUNCLEN_FWD  : ", TRUNCLEN_FWD)
log_r("  TRUNCLEN_REV  : ", TRUNCLEN_REV)
log_r("  MAX_EE_FWD    : ", MAX_EE_FWD)
log_r("  MAX_EE_REV    : ", MAX_EE_REV)
log_r("  TRUNC_Q       : ", TRUNC_Q)
log_r("  MAX_N         : ", MAX_N)

# Create temporary filtered directory
dir.create(FILTERED_DIR, recursive = TRUE, showWarnings = FALSE)

filt_fwd <- file.path(FILTERED_DIR, paste0(sample_names, "_R1_filt.fastq.gz"))
filt_rev <- file.path(FILTERED_DIR, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filt_fwd) <- sample_names
names(filt_rev) <- sample_names

out <- filterAndTrim(
    fwd        = fwd_files,
    filt       = filt_fwd,
    rev        = rev_files,
    filt.rev   = filt_rev,
    truncLen   = c(TRUNCLEN_FWD, TRUNCLEN_REV),
    maxEE      = c(MAX_EE_FWD,   MAX_EE_REV),
    truncQ     = TRUNC_Q,
    maxN       = MAX_N,
    compress   = COMPRESS,
    multithread = THREADS
)

rownames(out) <- sample_names
log_r("  filterAndTrim complete.")
log_r("  Reads in  : ", sum(out[, "reads.in"]))
log_r("  Reads out : ", sum(out[, "reads.out"]))

# --- Handle zero-read samples ------------------------------------------------
zero_samples <- sample_names[out[, "reads.out"] == 0]
if (length(zero_samples) > 0) {
    log_r("  WARNING: ", length(zero_samples), " sample(s) with zero reads after filtering:")
    for (s in zero_samples) log_r("    - ", s)
    pct_zero <- round(length(zero_samples) / length(sample_names) * 100, 1)
    if (pct_zero > 10) {
        log_r("  WARNING: ", pct_zero, "% of samples dropped — check TRUNCLEN and MAX_EE parameters.")
    }
}

# Keep only samples with reads
exists      <- out[, "reads.out"] > 0
filt_fwd    <- filt_fwd[exists]
filt_rev    <- filt_rev[exists]
sample_names_active <- sample_names[exists]

# --- Update read tracking ----------------------------------------------------
track_filter <- data.frame(
    step04_filtered = out[, "reads.out"],
    row.names       = sample_names
)
# Zero samples already have 0 from filterAndTrim
update_read_tracking(READ_TRACKING, track_filter)

# =============================================================================
# STEP 3 — LEARN ERROR RATES
# =============================================================================

log_r("----- Step 3: learnErrors : ", format(Sys.time(), "%H:%M:%S"), " -----")
log_r("  Threads (learnErrors) : ", THREADS_LEARN)

RDS_ERR_FWD <- file.path(STEP4_DIR, "err_fwd.rds")
RDS_ERR_REV <- file.path(STEP4_DIR, "err_rev.rds")

if (file.exists(RDS_ERR_FWD) && file.exists(RDS_ERR_REV)) {
    log_r("  Cached error models found — loading from RDS (delete to rerun learnErrors).")
    err_fwd <- readRDS(RDS_ERR_FWD)
    err_rev <- readRDS(RDS_ERR_REV)
} else {
    err_fwd <- learnErrors(filt_fwd, multithread = THREADS_LEARN)
    err_rev <- learnErrors(filt_rev, multithread = THREADS_LEARN)
    saveRDS(err_fwd, RDS_ERR_FWD)
    saveRDS(err_rev, RDS_ERR_REV)
    log_r("  Error models saved to RDS.")
}
log_r("  Error learning complete.")

# =============================================================================
# STEP 4 — DENOISE (DADA2)
# =============================================================================

log_r("----- Step 4: dada : ", format(Sys.time(), "%H:%M:%S"), " -----")
log_r("  Pooling : ", POOLING)

RDS_DADA_FWD <- file.path(STEP4_DIR, "dada_fwd.rds")
RDS_DADA_REV <- file.path(STEP4_DIR, "dada_rev.rds")

if (file.exists(RDS_DADA_FWD) && file.exists(RDS_DADA_REV)) {
    log_r("  Cached dada objects found — loading from RDS (delete to rerun dada).")
    dada_fwd <- readRDS(RDS_DADA_FWD)
    dada_rev <- readRDS(RDS_DADA_REV)
} else {
    dada_fwd <- dada(filt_fwd, err = err_fwd, pool = POOLING, multithread = THREADS)
    dada_rev <- dada(filt_rev, err = err_rev, pool = POOLING, multithread = THREADS)
    saveRDS(dada_fwd, RDS_DADA_FWD)
    saveRDS(dada_rev, RDS_DADA_REV)
    log_r("  Dada objects saved to RDS.")
}
log_r("  Denoising complete.")

# =============================================================================
# STEP 5 — MERGE PAIRED READS
# =============================================================================

log_r("----- Step 5: mergePairs : ", format(Sys.time(), "%H:%M:%S"), " -----")

# Dereplicate filtered reads — required for mergePairs regardless of pooling
derep_fwd <- derepFastq(filt_fwd)
derep_rev <- derepFastq(filt_rev)
names(derep_fwd) <- sample_names_active
names(derep_rev) <- sample_names_active

mergers <- mergePairs(
    dadaF     = dada_fwd,
    derepF    = derep_fwd,
    dadaR     = dada_rev,
    derepR    = derep_rev,
    minOverlap = MIN_OVERLAP,
    verbose   = FALSE
)

log_r("  Merging complete.")

# =============================================================================
# STEP 6 — SEQUENCE TABLE
# =============================================================================

log_r("----- Step 6: makeSequenceTable : ", format(Sys.time(), "%H:%M:%S"), " -----")

seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- sample_names_active
log_r("  ASVs before chimera removal : ", ncol(seqtab))
log_r("  Sequence length distribution :")
print(table(nchar(getSequences(seqtab))))

# =============================================================================
# STEP 7 — CHIMERA REMOVAL
# =============================================================================

log_r("----- Step 7: removeBimeraDenovo : ", format(Sys.time(), "%H:%M:%S"), " -----")
log_r("  Method : ", CHIMERA_METHOD)

seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method      = CHIMERA_METHOD,
    multithread = THREADS,
    verbose     = TRUE
)

log_r("  ASVs after chimera removal  : ", ncol(seqtab.nochim))
pct_retained <- round(sum(seqtab.nochim) / sum(seqtab) * 100, 1)
log_r("  Reads retained after chimera removal : ", pct_retained, "%")

# =============================================================================
# STEP 8 — SORT AND RENAME ASVs
# =============================================================================

log_r("----- Step 8: Sort and rename ASVs : ", format(Sys.time(), "%H:%M:%S"), " -----")

# Sort by decreasing total abundance, ties broken by sequence alphabetical order
abundance_order    <- order(colSums(seqtab.nochim),
                            colnames(seqtab.nochim),
                            decreasing = c(TRUE, FALSE),
                            method     = "radix")
seqtab.sorted      <- seqtab.nochim[, abundance_order, drop = FALSE]

# Generate zero-padded ASV IDs — number of digits determined by ASV count
n_asvs   <- ncol(seqtab.sorted)
n_digits <- nchar(as.character(n_asvs))
asv_ids  <- paste0("ASV_", formatC(seq(n_asvs), width = n_digits, flag = "0"))
colnames(seqtab.sorted) <- asv_ids

log_r("  Total ASVs : ", n_asvs)
log_r("  ASV ID format : ASV_", paste0(rep("0", n_digits - 1), collapse = ""), "1 ... ASV_", n_asvs)

# Restore zero-read samples as rows of zeros
if (length(zero_samples) > 0) {
    zero_mat <- matrix(0,
                       nrow     = length(zero_samples),
                       ncol     = n_asvs,
                       dimnames = list(zero_samples, asv_ids))
    seqtab.sorted <- rbind(seqtab.sorted, zero_mat)
    # Reorder to original sample order
    seqtab.sorted <- seqtab.sorted[sample_names, , drop = FALSE]
    log_r("  Zero-read samples restored as zero rows in community matrix.")
}

# =============================================================================
# STEP 9 — READ TRACKING UPDATE
# =============================================================================

log_r("----- Step 9: Update read tracking : ", format(Sys.time(), "%H:%M:%S"), " -----")

# Build per-sample tracking for active samples
get_n_reads <- function(dada_obj, names) {
    sapply(names, function(s) {
        idx <- which(names(dada_obj) == s)
        if (length(idx) == 0) return(0)
        sum(dada_obj[[idx]]$denoised)
    })
}

denoised_fwd <- sapply(dada_fwd, function(x) sum(x$denoised))
denoised_rev <- sapply(dada_rev, function(x) sum(x$denoised))
merged       <- sapply(mergers,   function(x) sum(x$abundance[x$accept]))
names(denoised_fwd) <- sample_names_active
names(denoised_rev) <- sample_names_active
names(merged)       <- sample_names_active

nochim <- rowSums(seqtab.sorted)

track_dada <- data.frame(
    step04_denoised_fwd = c(denoised_fwd, setNames(rep(0, length(zero_samples)), zero_samples))[sample_names],
    step04_denoised_rev = c(denoised_rev, setNames(rep(0, length(zero_samples)), zero_samples))[sample_names],
    step04_merged       = c(merged,       setNames(rep(0, length(zero_samples)), zero_samples))[sample_names],
    step04_nochim       = nochim[sample_names],
    row.names           = sample_names
)
update_read_tracking(READ_TRACKING, track_dada)
log_r("  Read tracking updated.")

# =============================================================================
# STEP 10 — WRITE OUTPUTS
# =============================================================================

log_r("----- Step 10: Write outputs : ", format(Sys.time(), "%H:%M:%S"), " -----")

# --- Community matrix --------------------------------------------------------
mr <- as.data.frame(seqtab.sorted)
write.table(mr,
            file      = OUT_MR,
            sep       = ";",
            quote     = FALSE,
            col.names = NA)
log_r("  Community matrix written : ", OUT_MR)

# --- FASTA -------------------------------------------------------------------
asv_sequences <- colnames(seqtab.nochim)[abundance_order]
write_fasta(
    sequences = asv_sequences,
    names     = asv_ids,
    file      = OUT_FASTA
)
log_r("  FASTA written : ", OUT_FASTA)

# --- seqtab RDS (sorted + renamed — reference for step 05) ------------------
saveRDS(seqtab.sorted, OUT_SEQTAB_RDS)
log_r("  seqtab RDS written : ", OUT_SEQTAB_RDS)

# =============================================================================
# STEP 11 — CLEANUP TEMPORARY FILES
# =============================================================================

log_r("----- Step 11: Cleanup : ", format(Sys.time(), "%H:%M:%S"), " -----")
unlink(FILTERED_DIR, recursive = TRUE)
log_r("  Temporary filtered directory removed.")

# =============================================================================
# FOOTER
# =============================================================================

log_r("=============================================================================")
log_r("  04_DADA2_ASV.R finished : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_r("  Log    : ", log_file)
log_r("  MR     : ", OUT_MR)
log_r("  FASTA  : ", OUT_FASTA)
log_r("  RDS    : ", OUT_SEQTAB_RDS)
log_r("=============================================================================")

sink()
