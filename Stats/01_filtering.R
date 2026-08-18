# =============================================================================
# File      : 01_filtering.R
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Author    : David Singer
# License   : MIT
# =============================================================================
#
# DESCRIPTION
# -----------
# Second script of the MarkerFlow statistics stage. Applies all matrix
# transformations to the raw data imported by 00_import.R: taxonomic filtering,
# rare ASV removal, quality filtering, sample exclusions, rarefaction, and
# residual correction. Produces a fully processed DS object saved back to
# <PROJECT>.RData. No figures are produced here — all diagnostics are in
# 02_dataset_quality.R.
#
# INPUT
# -----
# <PROJECT>.RData — produced by 00_import.R
#
# OUTPUT
# ------
# <PROJECT>.RData                            : updated — adds $MR_clean,
#                                              $ASS_clean, $MR_rare, $ASS_rare,
#                                              $MR_cor, $PARAMS per marker
# outputs/tables/<PROJECT>_01_journey_marker.csv : ASV/read counts at each
#                                              filtering step per marker
# outputs/tables/<PROJECT>_01_journey_sample.csv : reads and per-sample richness
#                                              at each filtering stage
# outputs/<PROJECT>_sessioninfo_01.txt       : R session information
#
# SLOTS ADDED TO DS PER MARKER
# --------------------------------
# $MR_clean   — after taxonomic + rare ASV + quality + sample exclusion filters
# $ASS_clean  — taxonomy table synced to MR_clean
# $MR_rare    — averaged rarefied MR_clean (RARE_ITER iterations, seed RARE_SEED)
# $ASS_rare   — taxonomy table synced to MR_rare
# $MR_cor     — named numeric vector: Observed ASV log-residuals (one per sample)
# $PARAMS     — all thresholds, decisions, and journey statistics
#
# $ASS (raw, from 00_import.R) is retained unchanged for full traceability.
# This stage is metadata-free — no sample metadata is used or required.
#
# SECTIONS
# --------
# 1. Parameters
# 2. Filtering functions
# 3. Apply filters per marker
# 4. Filtering journey tables
# 5. Rarefaction
# 6. Residual correction
# 7. Save DS.RData
# 8. Session info
# 9. Clean environment
#
# DEPENDENCIES
# ------------
# vegan >= 2.6
#
# REFERENCES
# ----------
# Callahan et al. (2016) DADA2. Nature Methods 13:581–583
# Oksanen et al. (2024) vegan. R package version 2.6-8
# Cameron & Tremblay (2020) Ecosphere 11:e03253
# =============================================================================


# =============================================================================
# SECTION 1 — PARAMETERS
# =============================================================================

# --- 1.1  Project name --------------------------------------------------------

PROJECT_NAME <- "MarkerFlow"

# --- 1.2  Paths ---------------------------------------------------------------

ROOT_DIR <- getwd()

PATH_RDATA       <- file.path(ROOT_DIR, paste0(PROJECT_NAME, ".RData"))
PATH_SESSION_TXT <- file.path(ROOT_DIR, "outputs",
                              paste0(PROJECT_NAME, "_sessioninfo_01.txt"))
PATH_JOURNEY_MARKER <- file.path(ROOT_DIR, "outputs", "tables",
                                 paste0(PROJECT_NAME,
                                        "_01_journey_marker.csv"))
PATH_JOURNEY_SAMPLE <- file.path(ROOT_DIR, "outputs", "tables",
                                 paste0(PROJECT_NAME,
                                        "_01_journey_sample.csv"))

# --- 1.3  Markers -------------------------------------------------------------
# Must match the markers imported in 00_import.R.

MARKERS <- c("18S", "ITS", "16S")

# --- 1.4  Taxonomic exclusion filter ------------------------------------------
# Per-marker named list of taxa to remove based on Division and Subdivision.
# Use list() for markers where no taxonomic exclusion is needed (the default).
# Applied only to markers where the list is non-empty.
#
# Example (18S / PR2 taxonomy) — to restrict an 18S dataset to protists by
# removing fungi, metazoa and plants, replace "18S" = list() with:
#   "18S" = list(
#     list(division = "Opisthokonta", subdivision = "Fungi"),
#     list(division = "Opisthokonta", subdivision = "Metazoa"),
#     list(division = "Streptophyta", subdivision = "Streptophyta_X")
#   )

EXCLUDE_TAXA <- list(
  "18S" = list(),   # no taxonomic exclusion by default (see example above)
  "ITS" = list(),
  "16S" = list()
)

# --- 1.5  Rare ASV filter -----------------------------------------------------
# An ASV is retained only if it has at least FILTER_MIN_READS reads in at
# least FILTER_MIN_SAMPLES samples. Removes sample-specific artefacts while
# preserving genuine rare taxa present across multiple samples.

FILTER_MIN_READS   <- 2   # minimum read count per ASV per sample
FILTER_MIN_SAMPLES <- 2   # minimum number of samples satisfying the above

# --- 1.6  Quality filter ------------------------------------------------------
# Removes ASVs with insufficient taxonomic assignment confidence.
# Thresholds apply to VSEARCH % identity and Domain-level bootstrap support.

FILTER_MIN_IDENTITY    <- 60    # minimum VSEARCH % identity to reference
FILTER_MIN_BOOT_DOMAIN <- 0.50  # minimum bootstrap support at Domain level

# --- 1.7  Sample exclusions ---------------------------------------------------
# Samples removed after visual inspection of read depth and diversity.
# Named list parallel to MARKERS. Use c() for markers with no exclusions.
# Document the reason for each excluded sample as an inline comment.

EXCLUDE_SAMPLES <- list(
  "18S" = c(),   # e.g. c("SampleA", "SampleB") — very low-read samples to drop
  "ITS" = c(),
  "16S" = c()
)

# --- 1.8  Rarefaction ---------------------------------------------------------
# Marker-specific rarefaction depth. Samples below threshold are excluded from
# MR_rare. Results are averaged over RARE_ITER iterations to reduce stochastic
# noise (Cameron & Tremblay 2020). Set from the read-depth diagnostics in
# 02_dataset_quality.R (and Sections 2.5 / 3.5 below).
# NOTE: the defaults below are LOW, sized for the small bundled example data.
# For a real dataset raise them (protist/fungal amplicons often 3,000–10,000).

RARE_THRESHOLDS <- list(
  "18S" = 1000,
  "ITS" = 1000,
  "16S" = 500
)
RARE_ITER <- 100   # number of rarefaction iterations to average
RARE_SEED <- 42    # random seed for reproducibility

# --- 1.9  Residual correction -------------------------------------------------
# Indices listed in INDICES_RESIDUAL are corrected for sequencing depth using
# linear regression residuals (index ~ f(reads)). All other indices are taken
# directly from MR_rare (rarefied values).
# Three transformations are tested (linear, log, sqrt); the best is selected
# by lowest mean |rho| between residuals and read depth across all indices.
# ALPHA_SIG controls whether a depth effect is considered significant enough
# to warrant correction.

INDICES_RESIDUAL <- "Observed"   # only Observed ASVs use residual correction
ALPHA_SIG        <- 0.05         # Spearman p-value threshold

# --- 1.10  ASS column name mapping --------------------------------------------
# Must match the values used in 00_import.R Section 1.6.

ASS_COL_ID       <- "ASV_ID"
ASS_COL_IDENTITY <- "Pct_identity"
ASS_COL_BOOT_DOM <- "Domain_boot"
ASS_COL_DIVISION <- "Division"
ASS_COL_SUBDIV   <- "Subdivision"

# --- 1.11  Marker display options ---------------------------------------------

MARKER_LABELS <- c(
  "18S" = "18S (protists)",
  "ITS" = "ITS (fungi)",
  "16S" = "16S (bacteria)"
)

MARKER_COLORS <- c(
  "18S" = "#2E86AB",
  "ITS" = "#E84855",
  "16S" = "#3BB273"
)

# --- 1.12  Figure output options ----------------------------------------------

SAVE_PDF  <- TRUE    # save figures as PDF
SAVE_JPEG <- TRUE    # save figures as JPEG (supplementary material)
FIG_DPI   <- 300     # JPEG resolution

PATH_FIG_DIR <- file.path(ROOT_DIR, "outputs", "figures")
dir.create(PATH_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1.13  Decision support parameters ---------------------------------------
# Used in Sections 2.5 and 3.5 only. Adjust after a first run.

N_LOW_SAMPLES <- 10     # number of lowest-read samples shown in bar plot
BS_QUANTILE   <- 0.15   # fraction of lowest-read samples used for broken
# stick analysis (0.15 = lower 15%; increase toward
# 0.50 if no break is detected)

# =============================================================================
# SECTION 2 — FUNCTIONS
# =============================================================================

# --- 2.1  Filtering functions --------------------------------------------------

#' Remove taxa listed in exclude_taxa based on Division and Subdivision.
#'
#' @param mr           data.frame. Samples x ASVs abundance matrix.
#' @param ass          data.frame. Taxonomic assignment table.
#' @param exclude_taxa List of lists, each with 'division' and 'subdivision'.
#' @param col_div      Character. Name of the Division column in ASS.
#' @param col_subdiv   Character. Name of the Subdivision column in ASS.
#' @return Filtered MR data.frame.
filter_taxa <- function(mr, ass, exclude_taxa,
                        col_div = "Division", col_subdiv = "Subdivision") {
  exclude <- rep(FALSE, nrow(ass))
  for (taxon in exclude_taxa) {
    exclude <- exclude | (
      !is.na(ass[[col_div]])    & ass[[col_div]]    == taxon$division &
        !is.na(ass[[col_subdiv]]) & ass[[col_subdiv]] == taxon$subdivision
    )
  }
  keep    <- rownames(ass)[!exclude]
  mr_filt <- mr[, colnames(mr) %in% keep, drop = FALSE]
  cat(sprintf("    Taxa filter     : %d -> %d ASVs | %d removed\n",
              ncol(mr), ncol(mr_filt), ncol(mr) - ncol(mr_filt)))
  for (taxon in exclude_taxa) {
    n <- sum(!is.na(ass[[col_div]])    & ass[[col_div]]    == taxon$division &
               !is.na(ass[[col_subdiv]]) & ass[[col_subdiv]] == taxon$subdivision)
    cat(sprintf("      - %s / %s : %d ASVs\n",
                taxon$division, taxon$subdivision, n))
  }
  return(mr_filt)
}

#' Remove rare ASVs below minimum abundance in minimum number of samples.
#'
#' An ASV is retained only if it has at least min_reads reads in at least
#' min_samples samples. Removes sample-specific sequencing artefacts while
#' preserving genuine rare taxa present across multiple samples.
#'
#' @param mr          data.frame. Samples x ASVs abundance matrix.
#' @param min_reads   Integer. Minimum read count per ASV per sample.
#' @param min_samples Integer. Minimum number of samples satisfying min_reads.
#' @return Filtered MR data.frame.
filter_rare_asvs <- function(mr, min_reads, min_samples) {
  keep       <- colSums(mr >= min_reads) >= min_samples
  reads_lost <- sum(mr) - sum(mr[, keep])
  cat(sprintf(
    "    Rare ASV filter : %d -> %d ASVs | %d removed (%.1f%%) | %d reads lost (%.1f%%)\n",
    ncol(mr), sum(keep), sum(!keep),
    sum(!keep) / ncol(mr) * 100,
    reads_lost, reads_lost / sum(mr) * 100))
  return(mr[, keep, drop = FALSE])
}

#' Remove ASVs below identity or domain bootstrap thresholds.
#'
#' ASVs with missing identity or bootstrap values are also removed, as these
#' indicate failed taxonomic assignment.
#'
#' @param mr              data.frame. Samples x ASVs abundance matrix.
#' @param ass             data.frame. Taxonomic assignment table.
#' @param min_identity    Numeric. Minimum VSEARCH % identity.
#' @param min_boot_domain Numeric. Minimum bootstrap at Domain level (0–1).
#' @param col_identity    Character. Name of the identity column in ASS.
#' @param col_boot_domain Character. Name of the Domain bootstrap column.
#' @return Filtered MR data.frame.
filter_quality <- function(mr, ass, min_identity, min_boot_domain,
                           col_identity    = "Pct_identity",
                           col_boot_domain = "Domain_boot") {
  asvs    <- colnames(mr)
  ass_sub <- ass[asvs, ]
  remove  <- (
    is.na(ass_sub[[col_identity]])    |
      ass_sub[[col_identity]]    < min_identity  |
      is.na(ass_sub[[col_boot_domain]]) |
      ass_sub[[col_boot_domain]] < min_boot_domain
  )
  reads_lost <- sum(mr[, remove])
  mr_filt    <- mr[, !remove, drop = FALSE]
  cat(sprintf(
    "    Quality filter  : %d -> %d ASVs | %d removed (%.1f%%) | %d reads lost (%.1f%%)\n",
    ncol(mr), ncol(mr_filt), sum(remove),
    sum(remove) / ncol(mr) * 100,
    reads_lost, reads_lost / sum(mr) * 100))
  cat(sprintf("      - Identity < %g%%     : %d ASVs\n",
              min_identity,
              sum(!is.na(ass_sub[[col_identity]]) &
                    ass_sub[[col_identity]] < min_identity)))
  cat(sprintf("      - Domain boot < %.2f : %d ASVs\n",
              min_boot_domain,
              sum(!is.na(ass_sub[[col_boot_domain]]) &
                    ass_sub[[col_boot_domain]] < min_boot_domain)))
  cat(sprintf("      - Missing values      : %d ASVs\n",
              sum(is.na(ass_sub[[col_identity]]) |
                    is.na(ass_sub[[col_boot_domain]]))))
  return(mr_filt)
}

#' Remove explicitly listed samples from an MR matrix.
#'
#' @param mr      data.frame. Samples x ASVs abundance matrix.
#' @param samples Character vector. Sample IDs to remove.
#' @param marker  Character. Marker name for warnings.
#' @return Filtered MR data.frame.
exclude_samples <- function(mr, samples, marker) {
  if (length(samples) == 0) {
    cat("    Sample exclusion: none\n")
    return(mr)
  }
  found   <- samples[samples %in% rownames(mr)]
  missing <- samples[!samples %in% rownames(mr)]
  if (length(missing) > 0)
    warning(marker, ": samples listed for exclusion but not found — ",
            paste(missing, collapse = "; "))
  mr_filt <- mr[!rownames(mr) %in% found, , drop = FALSE]
  cat(sprintf("    Sample exclusion: %d -> %d samples | removed: %s\n",
              nrow(mr), nrow(mr_filt), paste(found, collapse = "; ")))
  return(mr_filt)
}

#' Stop with an informative error if the filtered matrix has no rows or columns.
#'
#' @param mr     data.frame. Filtered MR matrix.
#' @param marker Character. Marker name for the error message.
#' @param step   Character. Filter step name for the error message.
check_not_empty <- function(mr, marker, step) {
  if (ncol(mr) == 0)
    stop(marker, ": no ASVs remaining after '", step,
         "'. Thresholds may be too stringent — check Section 1 parameters.")
  if (nrow(mr) == 0)
    stop(marker, ": no samples remaining after '", step,
         "'. Check EXCLUDE_SAMPLES and filter thresholds.")
}


# --- 2.2  Sync utilities ------------------------------------------------------

#' Subset ASS to ASVs present in the filtered MR, in the same column order.
#'
#' @param mr  data.frame. Filtered samples x ASVs matrix.
#' @param ass data.frame. Full taxonomic assignment table.
#' @return Subset of ass matching colnames(mr) in the same order.
sync_ass <- function(mr, ass) {
  ass[colnames(mr), , drop = FALSE]
}


# --- 2.3  Filtering journey tables --------------------------------------------

#' Append one filtering step to the marker-level filtering log.
#'
#' @param log       List. Current filtering log.
#' @param marker    Character. Marker name.
#' @param step      Character. Step label.
#' @param mr_before data.frame. MR immediately before this step.
#' @param mr_after  data.frame. MR immediately after this step.
#' @param note      Character. Free-text note (threshold; excluded IDs).
#' @return Updated filtering log list.
log_step <- function(log, marker, step, mr_before, mr_after, note = "") {
  log[[length(log) + 1]] <- list(
    marker       = marker,
    step         = step,
    samp_before  = nrow(mr_before),
    samp_after   = nrow(mr_after),
    asvs_before  = ncol(mr_before),
    asvs_after   = ncol(mr_after),
    reads_before = sum(mr_before),
    reads_after  = sum(mr_after),
    note         = note
  )
  return(log)
}

#' Build the marker-level journey table from the filtering log.
#'
#' Per-step percentages are relative to the state entering that step.
#' One FINAL row per marker compares raw import to the post-exclusion state.
#'
#' @param log List. Filtering log populated during Section 3.
#' @return A data.frame suitable for export as CSV.
build_journey_marker <- function(log) {
  
  rows <- lapply(log, function(x) {
    asvs_lost  <- x$asvs_before  - x$asvs_after
    reads_lost <- x$reads_before - x$reads_after
    data.frame(
      marker           = x$marker,
      step             = x$step,
      samples_before   = x$samp_before,
      samples_after    = x$samp_after,
      samples_removed  = x$samp_before  - x$samp_after,
      asvs_before      = x$asvs_before,
      asvs_after       = x$asvs_after,
      asvs_removed     = asvs_lost,
      asvs_removed_pct = round(asvs_lost  / x$asvs_before  * 100, 1),
      reads_before     = x$reads_before,
      reads_after      = x$reads_after,
      reads_lost       = reads_lost,
      reads_lost_pct   = round(reads_lost / x$reads_before * 100, 1),
      note             = x$note,
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, rows)
  
  # Append one FINAL row per marker: raw import vs post-exclusion state
  markers    <- unique(tbl$marker)
  final_rows <- lapply(markers, function(m) {
    raw  <- tbl[tbl$marker == m & tbl$step == "Raw", ]
    last <- tbl[tbl$marker == m, ]
    last <- last[nrow(last), ]
    asvs_lost  <- raw$asvs_before  - last$asvs_after
    reads_lost <- raw$reads_before - last$reads_after
    data.frame(
      marker           = m,
      step             = "FINAL",
      samples_before   = raw$samples_before,
      samples_after    = last$samples_after,
      samples_removed  = raw$samples_before  - last$samples_after,
      asvs_before      = raw$asvs_before,
      asvs_after       = last$asvs_after,
      asvs_removed     = asvs_lost,
      asvs_removed_pct = round(asvs_lost  / raw$asvs_before  * 100, 1),
      reads_before     = raw$reads_before,
      reads_after      = last$reads_after,
      reads_lost       = reads_lost,
      reads_lost_pct   = round(reads_lost / raw$reads_before * 100, 1),
      note             = "raw vs post-exclusion",
      stringsAsFactors = FALSE
    )
  })
  
  tbl_final <- do.call(rbind, lapply(markers, function(m) {
    rbind(tbl[tbl$marker == m, ], final_rows[[match(m, markers)]])
  }))
  rownames(tbl_final) <- NULL
  return(tbl_final)
}

#' Build the per-sample journey table across filtering stages.
#'
#' Records reads and observed ASV richness (non-zero columns per sample row)
#' at five stages: raw, taxa filter, rare ASV filter, quality filter, and
#' rarefaction. Samples excluded before rarefaction are retained in the table
#' with "excluded" in the rarefied columns for full traceability.
#' For markers with no taxonomic filter (ITS, 16S), taxa columns equal raw.
#'
#' @param mr_raw    data.frame. Raw MR from 00_import.R.
#' @param mr_taxa   data.frame. MR after taxonomic filter (or mr_raw if none).
#' @param mr_rareASV data.frame. MR after rare ASV filter.
#' @param mr_qual   data.frame. MR after quality filter (= MR_clean before
#'                  sample exclusions, used for per-sample per-step tracking).
#' @param mr_rare   data.frame. Rarefied MR (MR_rare — may have fewer samples).
#' @return A data.frame with one row per sample (all samples from mr_raw).

build_journey_sample <- function(mr_raw, mr_taxa, mr_rareASV,
                                 mr_qual, mr_rare) {
  
  # Vectorised stats for one matrix: reads and richness per sample row
  stage_stats <- function(mr, all_samples) {
    reads <- rowSums(mr)
    asvs  <- rowSums(mr > 0)
    data.frame(
      reads = reads[all_samples],
      asvs  = asvs[all_samples],
      row.names = all_samples
    )
  }
  
  all_samples <- rownames(mr_raw)
  
  raw     <- stage_stats(mr_raw,     all_samples)
  taxa    <- stage_stats(mr_taxa,    all_samples)
  rareASV <- stage_stats(mr_rareASV, all_samples)
  qual    <- stage_stats(mr_qual,    all_samples)
  
  # Rarefied stage — only samples present in mr_rare get values
  rare_reads <- rare_asvs <- rep(NA_integer_, length(all_samples))
  names(rare_reads) <- names(rare_asvs) <- all_samples
  in_rare <- intersect(all_samples, rownames(mr_rare))
  rare_reads[in_rare] <- as.integer(rowSums(mr_rare)[in_rare])
  rare_asvs[in_rare]  <- as.integer(rowSums(mr_rare > 0)[in_rare])
  
  data.frame(
    sample         = all_samples,
    raw_reads      = raw$reads,
    raw_asvs       = raw$asvs,
    taxa_reads     = taxa$reads,
    taxa_asvs      = taxa$asvs,
    rareASV_reads  = rareASV$reads,
    rareASV_asvs   = rareASV$asvs,
    quality_reads  = qual$reads,
    quality_asvs   = qual$asvs,
    rarefied_reads = ifelse(is.na(rare_reads), "excluded",
                            as.character(rare_reads)),
    rarefied_asvs  = ifelse(is.na(rare_asvs),  "excluded",
                            as.character(rare_asvs)),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

# --- 2.4  Rarefaction ---------------------------------------------------------

#' Rarefy a matrix to even depth using averaged multiple iterations.
#'
#' Samples below threshold are excluded before rarefaction. Each sample is
#' rarefied iter times and results are averaged to reduce stochastic noise
#' from single subsampling (Cameron & Tremblay 2020). ASVs with zero counts
#' after averaging are removed.
#'
#' @param mr        data.frame. Samples x ASVs matrix (MR_clean).
#' @param threshold Integer. Rarefaction depth.
#' @param iter      Integer. Number of rarefaction iterations to average.
#' @param seed      Integer. Random seed for reproducibility.
#' @param marker    Character. Marker name for console output.
#' @return A rarefied data.frame (samples x ASVs), averaged over iterations.
rarefy_average <- function(mr, threshold, iter, seed, marker) {
  
  reads   <- rowSums(mr)
  n_below <- sum(reads < threshold)
  
  if (n_below > 0) {
    cat(sprintf(
      "  %-6s | %d sample(s) below threshold (%s reads) — excluded from MR_rare:\n",
      marker, n_below, format(threshold, big.mark = ",")))
    for (s in names(reads[reads < threshold]))
      cat(sprintf("           |   %s : %s reads\n",
                  s, format(reads[s], big.mark = ",")))
    mr <- mr[reads >= threshold, , drop = FALSE]
  }
  
  cat(sprintf(
    "  %-6s | Rarefying %d samples to %s reads (%d iterations, seed = %d) ...\n",
    marker, nrow(mr), format(threshold, big.mark = ","), iter, seed))
  
  set.seed(seed)
  rare_sum <- matrix(0,
                     nrow     = nrow(mr),
                     ncol     = ncol(mr),
                     dimnames = dimnames(mr))
  for (i in seq_len(iter))
    rare_sum <- rare_sum + vegan::rrarefy(mr, sample = threshold)
  
  # Average and round to nearest integer.
  # Note: rounding introduces minor deviations from the exact threshold depth
  # (~1 read per sample); Shannon and Simpson are robust to this deviation.
  mr_rare              <- as.data.frame(round(rare_sum / iter))
  rownames(mr_rare)    <- rownames(mr)  
  # Remove ASVs with zero counts after averaging
  mr_rare <- mr_rare[, colSums(mr_rare) > 0, drop = FALSE]
  
  cat(sprintf("  %-6s | MR_rare : %d samples x %d ASVs\n",
              marker, nrow(mr_rare), ncol(mr_rare)))
  return(mr_rare)
}


# --- 2.5  Residual correction -------------------------------------------------

#' Test Spearman correlation between diversity index values and read depth.
#'
#' @param values Numeric vector. Diversity index values.
#' @param reads  Numeric vector. Total reads per sample.
#' @return Named numeric vector: rho and p.
spearman_depth <- function(values, reads) {
  test <- cor.test(values, reads, method = "spearman", exact = FALSE)
  c(rho = unname(round(test$estimate, 3)),
    p   = unname(round(test$p.value,  4)))
}

#' Compute residuals of index ~ f(reads) for one transformation.
#'
#' Uses na.exclude to preserve output vector length when reads contains NA.
#'
#' @param values    Numeric vector. Diversity index values.
#' @param reads     Numeric vector. Total reads per sample.
#' @param transform Character. One of "linear", "log", "sqrt".
#' @return Numeric vector of residuals, same length as values.
compute_residuals <- function(values, reads, transform) {
  x <- switch(transform,
              linear = reads,
              log    = log(reads),
              sqrt   = sqrt(reads))
  resids <- residuals(lm(values ~ x, na.action = na.exclude))
  out    <- rep(NA_real_, length(values))
  out[seq_along(resids)[!is.na(resids)]] <- resids[!is.na(resids)]
  return(out)
}

#' Evaluate all transformations for one index x marker.
#'
#' @param values  Numeric vector. Index values from MR_clean (non-rarefied).
#' @param reads   Numeric vector. Total reads per sample from MR_clean.
#' @param marker  Character. Marker name.
#' @param index   Character. Index name.
#' @param alpha   Numeric. Significance threshold for depth effect.
#' @return A data.frame with one row per transformation.
evaluate_transforms <- function(values, reads, marker, index, alpha = 0.05) {
  before <- spearman_depth(values, reads)
  rows <- lapply(c("linear", "log", "sqrt"), function(tr) {
    resids <- compute_residuals(values, reads, tr)
    after  <- spearman_depth(resids, reads)
    data.frame(
      marker             = marker,
      index              = index,
      transform          = tr,
      rho_before         = before["rho"],
      p_before           = before["p"],
      rho_after          = after["rho"],
      p_after            = after["p"],
      correction_applied = before["p"] <= alpha,
      validated          = after["p"]  >  alpha,
      stringsAsFactors   = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Select the best transformation: lowest mean |rho_after| across all indices.
#'
#' @param eval_tbl data.frame. Output of evaluate_transforms() for all indices.
#' @return Character. Name of the best transformation.
select_best_transform <- function(eval_tbl) {
  avg_rho <- tapply(abs(eval_tbl$rho_after), eval_tbl$transform,
                    mean, na.rm = TRUE)
  names(which.min(avg_rho))
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("--- Loading DS.RData ------------------------------\n")
if (!file.exists(PATH_RDATA))
  stop("DS.RData not found at: ", PATH_RDATA,
       "\nRun 00_import.R first.")
load(PATH_RDATA)
cat("  Loaded:", PATH_RDATA, "\n")
cat("  Markers found:", paste(names(DS), collapse = "; "), "\n\n")

# =============================================================================
# SECTION 2.5 — DECISION SUPPORT: PRE-FILTERING DIAGNOSTICS (optional)
# =============================================================================
#
# Runs on raw matrices (DS[[marker]]$MR) before any filtering.
# Use these figures to guide EXCLUDE_SAMPLES and RARE_THRESHOLDS in Section 1,
# then re-run from Section 3.
#
# Figures produced (per marker, combined into one multi-panel PDF per figure):
#   DS_01_25_lowest_reads  : bar plot of N_LOW_SAMPLES lowest-read samples
#   DS_01_25_reads_vs_asvs : raw reads vs observed ASVs (Spearman rho)
#   DS_01_25_rare_support  : sorted reads dot plot + broken stick suggestion
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 2.5: Pre-filtering decision support\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Figures saved to:", PATH_FIG_DIR, "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
for (pkg in c("ggplot2", "patchwork"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(ggplot2)
library(patchwork)

# --- save_figure() helper -----------------------------------------------------
# Shared across Sections 2.5, 3.5 and 5.5.

#' Save a ggplot object as PDF and/or JPEG.
#'
#' @param fig      ggplot. Figure to save.
#' @param name     Character. File name stem (no extension).
#' @param width    Numeric. Width in inches.
#' @param height   Numeric. Height in inches.
save_figure <- function(fig, name, width = 12, height = 6) {
  base <- file.path(PATH_FIG_DIR, paste0(PROJECT_NAME, "_", name))
  if (SAVE_PDF)
    ggplot2::ggsave(paste0(base, ".pdf"),  plot = fig,
                    width = width, height = height)
  if (SAVE_JPEG)
    ggplot2::ggsave(paste0(base, ".jpg"),  plot = fig,
                    width = width, height = height, dpi = FIG_DPI)
  cat(sprintf("  Saved: %s [PDF=%s JPEG=%s]\n",
              basename(base), SAVE_PDF, SAVE_JPEG))
}

# =============================================================================
# 2.5.1 — LOWEST-READ SAMPLES BAR PLOT
# =============================================================================
# Identifies samples with very low read counts that are candidates for manual
# exclusion via EXCLUDE_SAMPLES in Section 1.
# =============================================================================

cat("--- 2.5.1  Lowest-read samples -----------------------\n\n")

fig_low <- setNames(lapply(MARKERS, function(marker) {
  
  mr     <- DS[[marker]]$MR
  reads  <- sort(rowSums(mr))[seq_len(min(N_LOW_SAMPLES, nrow(mr)))]
  df     <- data.frame(
    sample = factor(names(reads), levels = names(reads)),
    reads  = as.integer(reads)
  )
  
  cat(sprintf("  %s — %d lowest-read samples:\n", marker, nrow(df)))
  for (i in seq_len(nrow(df)))
    cat(sprintf("    %-12s : %s reads\n",
                df$sample[i], format(df$reads[i], big.mark = ",")))
  cat("\n")
  
  ggplot(df, aes(x = sample, y = reads)) +
    geom_col(fill = MARKER_COLORS[marker], alpha = 0.85) +
    geom_text(aes(label = format(reads, big.mark = ",")),
              hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                       labels = function(x) format(x, big.mark = ",")) +
    labs(title    = MARKER_LABELS[marker],
         subtitle = paste0("Lowest ", N_LOW_SAMPLES, " samples by read count"),
         x        = NULL,
         y        = "Total reads") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          panel.grid.major.y = element_blank())
  
}), MARKERS)

fig_low_combined <- wrap_plots(fig_low, ncol = length(MARKERS)) +
  plot_annotation(
    title    = paste0(PROJECT_NAME, " — Lowest-read samples (raw matrix)"),
    subtitle = paste0("Guide for EXCLUDE_SAMPLES in Section 1 | ",
                      "N = ", N_LOW_SAMPLES, " per marker"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

print(fig_low_combined)
save_figure(fig_low_combined, "01_25_lowest_reads",
            width = 5 * length(MARKERS), height = 6)

# =============================================================================
# 2.5.2 — RAW READS vs OBSERVED ASVs
# =============================================================================
# Assesses depth-richness correlation on the raw matrix.
# A significant Spearman rho suggests residual correction will be needed.
# =============================================================================

cat("--- 2.5.2  Raw reads vs observed ASVs ----------------\n\n")

fig_corr <- setNames(lapply(MARKERS, function(marker) {
  
  mr      <- DS[[marker]]$MR
  reads   <- rowSums(mr)
  asvs    <- vegan::specnumber(round(mr))
  test    <- cor.test(reads, asvs, method = "spearman", exact = FALSE)
  rho     <- round(test$estimate, 3)
  p_val   <- test$p.value
  p_label <- ifelse(p_val < 0.001, "p < 0.001",
                    paste0("p = ", round(p_val, 3)))
  sig     <- p_val <= ALPHA_SIG
  
  cat(sprintf("  %s : Spearman rho = %s | %s | %s\n",
              marker, rho, p_label,
              ifelse(sig, "significant — correction likely needed",
                     "not significant")))
  
  df <- data.frame(reads = reads, asvs = asvs)
  
  ggplot(df, aes(x = reads, y = asvs)) +
    geom_point(color = MARKER_COLORS[marker], alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE,
                color = ifelse(sig, "#CC0000", "grey40"),
                linewidth = 0.8) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(title    = MARKER_LABELS[marker],
         subtitle = paste0("Spearman rho = ", rho, " | ", p_label),
         x        = "Raw reads per sample",
         y        = "Observed ASVs") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(
            color = ifelse(sig, "#CC0000", "grey40")))
  
}), MARKERS)

cat("\n")

fig_corr_combined <- wrap_plots(fig_corr, ncol = length(MARKERS)) +
  plot_annotation(
    title    = paste0(PROJECT_NAME,
                      " — Raw reads vs observed ASVs (raw matrix)"),
    subtitle = paste0("Red regression line = significant depth effect ",
                      "(p <= ", ALPHA_SIG, ") | ",
                      "Guide for residual correction decision"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

print(fig_corr_combined)
save_figure(fig_corr_combined, "01_25_reads_vs_asvs",
            width = 5 * length(MARKERS), height = 5)

# =============================================================================
# 2.5.3 — SORTED READ DEPTH + BROKEN STICK
# =============================================================================
# Sorted dot plot of read depth per sample with the current RARE_THRESHOLDS
# overlaid. Broken stick analysis on the lower BS_QUANTILE fraction suggests
# a natural break point for the rarefaction threshold decision.
# =============================================================================

cat("--- 2.5.3  Sorted read depth + broken stick ----------\n\n")

#' Fit a broken stick model to find the optimal breakpoint in sorted reads.
#'
#' Evaluates all possible breakpoints within the lower quantile of the read
#' depth distribution. The breakpoint minimising the total residual sum of
#' squares from two fitted line segments is returned as the suggested
#' rarefaction threshold.
#'
#' @param reads     Numeric vector. Per-sample read counts.
#' @param quantile  Numeric. Upper quantile limit for breakpoint search (0–1).
#' @return A list: $threshold (suggested depth), $break_idx (position in
#'         sorted reads), $rss (residual sum of squares at best break).
broken_stick <- function(reads, quantile = 0.15) {
  
  sorted  <- sort(reads)
  n       <- length(sorted)
  n_lower <- max(4, floor(n * quantile))
  x       <- seq_len(n)
  
  best_rss <- Inf
  best_idx <- NA
  
  for (bp in 2:(n_lower - 1)) {
    left  <- x[1:bp]
    right <- x[bp:n]
    rss   <- sum(residuals(lm(sorted[1:bp]   ~ left))^2) +
      sum(residuals(lm(sorted[bp:n]   ~ right))^2)
    if (rss < best_rss) {
      best_rss <- rss
      best_idx <- bp
    }
  }
  
  list(threshold = sorted[best_idx],
       break_idx = best_idx,
       rss       = best_rss)
}

fig_rare <- setNames(lapply(MARKERS, function(marker) {
  
  mr      <- DS[[marker]]$MR
  reads   <- sort(rowSums(mr))
  n       <- length(reads)
  df      <- data.frame(rank  = seq_len(n),
                        reads = reads,
                        label = names(reads))
  
  bs        <- broken_stick(reads, quantile = BS_QUANTILE)
  threshold <- RARE_THRESHOLDS[[marker]]
  
  cat(sprintf("  %s:\n", marker))
  cat(sprintf("    Current threshold    : %s reads\n",
              format(threshold, big.mark = ",")))
  cat(sprintf("    Broken stick suggest.: %s reads (rank %d of %d)\n\n",
              format(bs$threshold, big.mark = ","), bs$break_idx, n))
  
  # Samples below current threshold
  n_below <- sum(reads < threshold)
  
  ggplot(df, aes(x = rank, y = reads)) +
    geom_point(
      aes(color = reads < threshold),
      size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = MARKER_COLORS[marker], "TRUE" = "#CC0000"),
      labels = c("FALSE" = "Retained", "TRUE"  = "Below threshold"),
      name   = NULL) +
    geom_hline(yintercept = threshold,
               linetype = "dashed", color = "#CC0000", linewidth = 0.8) +
    geom_hline(yintercept = bs$threshold,
               linetype = "dotted", color = "grey40", linewidth = 0.8) +
    annotate("text",
             x = n * 0.98, y = threshold,
             label = paste0("Threshold: ",
                            format(threshold, big.mark = ",")),
             hjust = 1, vjust = -0.5,
             color = "#CC0000", size = 3) +
    annotate("text",
             x = n * 0.98, y = bs$threshold,
             label = paste0("Broken stick: ",
                            format(bs$threshold, big.mark = ",")),
             hjust = 1, vjust = -0.5,
             color = "grey40", size = 3) +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(title    = MARKER_LABELS[marker],
         subtitle = paste0(n_below, " sample(s) below threshold | ",
                           "Broken stick quantile = ", BS_QUANTILE),
         x        = "Sample rank (sorted by read count)",
         y        = "Total reads") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold"),
          legend.position = "bottom")
  
}), MARKERS)

fig_rare_combined <- wrap_plots(fig_rare, ncol = length(MARKERS)) +
  plot_annotation(
    title    = paste0(PROJECT_NAME,
                      " — Sorted read depth (raw matrix)"),
    subtitle = paste0(
      "Red dashed = current RARE_THRESHOLDS | ",
      "Grey dotted = broken stick suggestion | ",
      "Red points = samples excluded from rarefaction"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

print(fig_rare_combined)
save_figure(fig_rare_combined, "01_25_rare_support",
            width = 5 * length(MARKERS), height = 6)

cat("--- Section 2.5 complete -----------------------------\n")
cat("  Review figures in:", PATH_FIG_DIR, "\n")
cat("  Adjust EXCLUDE_SAMPLES and RARE_THRESHOLDS in Section 1 if needed,\n")
cat("  then re-run from Section 3.\n\n")


# Review the figures above, then adjust EXCLUDE_SAMPLES / RARE_THRESHOLDS in
# Section 1 if needed and re-run from Section 3. (Parameters live only in
# Section 1 — no mid-script overrides.)



# =============================================================================
# SECTION 3 — APPLY FILTERS PER MARKER
# =============================================================================

# Load DS.RData produced by 00_import.R
cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Filtering\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Project root :", ROOT_DIR, "\n")
cat("  Markers      :", paste(MARKERS, collapse = "; "), "\n")
cat("=============================================================================\n\n")

cat("--- Loading DS.RData ------------------------------\n")
if (!file.exists(PATH_RDATA))
  stop("DS.RData not found at: ", PATH_RDATA,
       "\nRun 00_import.R first.")
load(PATH_RDATA)
cat("  Loaded:", PATH_RDATA, "\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
library(vegan)

# Initialise filtering log and per-stage matrix storage
# STAGE_MR stores intermediate matrices per marker for the per-sample
# journey table (Section 4). Keys: raw, taxa, rareASV, quality.
# The rarefied stage is added in Section 5 after rarefy_average().
FILTER_LOG <- list()
STAGE_MR   <- list()

# --- Loop over markers --------------------------------------------------------
for (marker in MARKERS) {
  
  cat(sprintf("--- %s -------------------------------------------\n", marker))
  
  mr  <- DS[[marker]]$MR
  ass <- DS[[marker]]$ASS
  
  # Record raw state
  STAGE_MR[[marker]] <- list(raw = mr)
  FILTER_LOG <- log_step(FILTER_LOG, marker, "Raw",
                         mr_before = mr, mr_after = mr)
  
  # --- 3.1  Taxonomic filter (18S only) ---------------------------------------
  if (length(EXCLUDE_TAXA[[marker]]) > 0) {
    mr_prev <- mr
    mr      <- filter_taxa(mr_prev, ass,
                           exclude_taxa = EXCLUDE_TAXA[[marker]],
                           col_div      = ASS_COL_DIVISION,
                           col_subdiv   = ASS_COL_SUBDIV)
    check_not_empty(mr, marker, "Taxa filter")
    FILTER_LOG <- log_step(FILTER_LOG, marker, "Taxa filter",
                           mr_before = mr_prev,
                           mr_after  = mr,
                           note      = paste(
                             sapply(EXCLUDE_TAXA[[marker]],
                                    function(t) t$subdivision),
                             collapse = "; "))
  }
  STAGE_MR[[marker]]$taxa <- mr   # equals raw for ITS and 16S
  
  # --- 3.2  Rare ASV filter ---------------------------------------------------
  mr_prev <- mr
  mr      <- filter_rare_asvs(mr_prev, FILTER_MIN_READS, FILTER_MIN_SAMPLES)
  check_not_empty(mr, marker, "Rare ASV filter")
  FILTER_LOG <- log_step(FILTER_LOG, marker, "Rare ASV filter",
                         mr_before = mr_prev,
                         mr_after  = mr,
                         note      = paste0("min ", FILTER_MIN_READS,
                                            " reads in min ",
                                            FILTER_MIN_SAMPLES, " samples"))
  STAGE_MR[[marker]]$rareASV <- mr
  
  # --- 3.3  Quality filter ----------------------------------------------------
  mr_prev <- mr
  mr      <- filter_quality(mr_prev, ass,
                            min_identity    = FILTER_MIN_IDENTITY,
                            min_boot_domain = FILTER_MIN_BOOT_DOMAIN,
                            col_identity    = ASS_COL_IDENTITY,
                            col_boot_domain = ASS_COL_BOOT_DOM)
  check_not_empty(mr, marker, "Quality filter")
  FILTER_LOG <- log_step(FILTER_LOG, marker, "Quality filter",
                         mr_before = mr_prev,
                         mr_after  = mr,
                         note      = paste0("identity >= ",
                                            FILTER_MIN_IDENTITY,
                                            "%; domain boot >= ",
                                            FILTER_MIN_BOOT_DOMAIN))
  STAGE_MR[[marker]]$quality <- mr
  
  # --- 3.4  Sample exclusions -------------------------------------------------
  mr_prev <- mr
  mr      <- exclude_samples(mr_prev, EXCLUDE_SAMPLES[[marker]], marker)
  check_not_empty(mr, marker, "Sample exclusion")
  FILTER_LOG <- log_step(FILTER_LOG, marker, "Sample exclusion",
                         mr_before = mr_prev,
                         mr_after  = mr,
                         note      = ifelse(
                           length(EXCLUDE_SAMPLES[[marker]]) > 0,
                           paste(EXCLUDE_SAMPLES[[marker]], collapse = "; "),
                           "none"))
  
  # --- Sync ASS_clean to MR_clean ---------------------------------------------
  ass_clean <- sync_ass(mr, ass)
  
  # --- Store MR_clean and ASS_clean in DS ----------------------------------
  DS[[marker]]$MR_clean  <- mr
  DS[[marker]]$ASS_clean <- ass_clean
  cat(sprintf("  MR_clean  : %d samples x %d ASVs | %s reads\n",
              nrow(mr), ncol(mr), format(sum(mr), big.mark = ",")))
  cat(sprintf("  ASS_clean : %d ASVs\n\n", nrow(ass_clean)))
  
}  # end marker loop

cat("--- Marker-level journey table -----------------------\n\n")

JOURNEY_MARKER <- build_journey_marker(FILTER_LOG)

write.table(JOURNEY_MARKER,
            file      = PATH_JOURNEY_MARKER,
            sep       = ";",
            row.names = FALSE,
            quote     = TRUE)

cat(sprintf("  %-6s | %-20s | %8s | %8s | %9s | %12s\n",
            "Marker", "Step", "Samples", "ASVs", "ASVs %", "Reads"))
cat(sprintf("  %s\n", strrep("-", 72)))

for (i in seq_len(nrow(JOURNEY_MARKER))) {
  r <- JOURNEY_MARKER[i, ]
  cat(sprintf("  %-6s | %-20s | %8d | %8d | %8.1f%% | %12s\n",
              r$marker, r$step,
              r$samples_after,
              r$asvs_after,
              100 - r$asvs_removed_pct,
              format(r$reads_after, big.mark = ",")))
}

cat(sprintf("\n  Saved to: %s\n\n", PATH_JOURNEY_MARKER))

cat("  Per-sample journey table: built in Section 5 after rarefaction.\n\n")

# =============================================================================
# SECTION 3.5 — DECISION SUPPORT: POST-FILTERING RAREFACTION SUPPORT (optional)
# =============================================================================
#
# Runs on MR_clean — after all filters and sample exclusions.
# This is the definitive read depth distribution for the rarefaction threshold
# decision. Produces one figure per marker combined into a single PDF:
#   DS_01_35_rare_support : sorted reads dot plot + broken stick suggestion
#
# Review the figure, adjust RARE_THRESHOLDS below if needed, then re-run
# from Section 4.
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 3.5: Post-filtering rarefaction support\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

cat("--- Sorted read depth (MR_clean) ---------------------\n\n")

fig_rare_clean <- setNames(lapply(MARKERS, function(marker) {
  
  mr    <- DS[[marker]]$MR_clean
  reads <- sort(rowSums(mr))
  n     <- length(reads)
  df    <- data.frame(rank  = seq_len(n),
                      reads = reads,
                      label = names(reads))
  
  bs        <- broken_stick(reads, quantile = BS_QUANTILE)
  threshold <- RARE_THRESHOLDS[[marker]]
  n_below   <- sum(reads < threshold)
  
  cat(sprintf("  %s:\n", marker))
  cat(sprintf("    Samples in MR_clean  : %d\n", n))
  cat(sprintf("    Current threshold    : %s reads\n",
              format(threshold, big.mark = ",")))
  cat(sprintf("    Broken stick suggest.: %s reads (rank %d of %d)\n",
              format(bs$threshold, big.mark = ","), bs$break_idx, n))
  cat(sprintf("    Samples below threshold: %d\n\n", n_below))
  # Exact read counts around the threshold
  n_show  <- 10
  below   <- tail(reads[reads <  threshold], n_show)
  above   <- head(reads[reads >= threshold], n_show)
  nearest <- c(below, above)
  
  cat(sprintf("    Samples nearest to threshold (%s reads):\n",
              format(threshold, big.mark = ",")))
  for (s in names(nearest))
    cat(sprintf("      %-12s : %s reads  [%s]\n",
                s,
                format(nearest[s], big.mark = ","),
                ifelse(nearest[s] < threshold, "below", "above")))
  cat("\n")
  ggplot(df, aes(x = rank, y = reads)) +
    geom_point(
      aes(color = reads < threshold),
      size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("FALSE" = MARKER_COLORS[marker], "TRUE" = "#CC0000"),
      labels = c("FALSE" = "Retained", "TRUE"  = "Below threshold"),
      name   = NULL) +
    geom_hline(yintercept = threshold,
               linetype = "dashed", color = "#CC0000", linewidth = 0.8) +
    geom_hline(yintercept = bs$threshold,
               linetype = "dotted", color = "grey40", linewidth = 0.8) +
    annotate("text",
             x = n * 0.98, y = threshold,
             label = paste0("Threshold: ",
                            format(threshold, big.mark = ",")),
             hjust = 1, vjust = -0.5,
             color = "#CC0000", size = 3) +
    annotate("text",
             x = n * 0.98, y = bs$threshold,
             label = paste0("Broken stick: ",
                            format(bs$threshold, big.mark = ",")),
             hjust = 1, vjust = -0.5,
             color = "grey40", size = 3) +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(title    = MARKER_LABELS[marker],
         subtitle = paste0(n_below, " sample(s) below threshold | ",
                           "Broken stick quantile = ", BS_QUANTILE),
         x        = "Sample rank (sorted by read count)",
         y        = "Total reads (MR_clean)") +
    theme_bw(base_size = 10) +
    theme(plot.title      = element_text(face = "bold"),
          legend.position = "bottom")
  
}), MARKERS)

fig_rare_clean_combined <- wrap_plots(fig_rare_clean, ncol = length(MARKERS)) +
  plot_annotation(
    title    = paste0(PROJECT_NAME,
                      " — Sorted read depth (MR_clean — post-filtering)"),
    subtitle = paste0(
      "Red dashed = current RARE_THRESHOLDS | ",
      "Grey dotted = broken stick suggestion | ",
      "Red points = samples excluded from rarefaction"),
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

print(fig_rare_clean_combined)
save_figure(fig_rare_clean_combined, "01_35_rare_support",
            width = 5 * length(MARKERS), height = 6)

cat("--- Section 3.5 complete -----------------------------\n")
cat("  Review figure in:", PATH_FIG_DIR, "\n\n")

# Review the figure above, then adjust RARE_THRESHOLDS in Section 1 if needed
# and re-run from Section 4. (Parameters live only in Section 1.)

# =============================================================================
# SECTION 4 — RAREFACTION
# =============================================================================

cat("--- Rarefaction --------------------------------------\n\n")

for (marker in MARKERS) {
  
  mr_rare <- rarefy_average(
    mr        = DS[[marker]]$MR_clean,
    threshold = RARE_THRESHOLDS[[marker]],
    iter      = RARE_ITER,
    seed      = RARE_SEED,
    marker    = marker
  )
  
  DS[[marker]]$MR_rare <- mr_rare
  # Sync ASS_clean to MR_rare — rarefaction may remove zero-count ASV columns
  DS[[marker]]$ASS_rare <- DS[[marker]]$ASS_clean[
    rownames(DS[[marker]]$ASS_clean) %in% colnames(mr_rare), , drop = FALSE]
  
  cat(sprintf("  ASS_rare : %d ASVs (vs %d in ASS_clean — %d removed by rarefaction)\n",
              nrow(DS[[marker]]$ASS_rare),
              nrow(DS[[marker]]$ASS_clean),
              nrow(DS[[marker]]$ASS_clean) - nrow(DS[[marker]]$ASS_rare)))
}  # end marker loop

# --- 4.1  Per-sample journey table --------------------------------------------
# Built here once MR_rare is available for all markers.
# Samples excluded during rarefaction appear as "excluded" in rarefied columns.

cat("\n--- Per-sample journey table -------------------------\n\n")

JOURNEY_SAMPLE <- do.call(rbind, lapply(MARKERS, function(marker) {
  tbl <- build_journey_sample(
    mr_raw     = STAGE_MR[[marker]]$raw,
    mr_taxa    = STAGE_MR[[marker]]$taxa,
    mr_rareASV = STAGE_MR[[marker]]$rareASV,
    mr_qual    = STAGE_MR[[marker]]$quality,
    mr_rare    = DS[[marker]]$MR_rare
  )
  tbl$marker <- marker
  tbl        <- tbl[, c("marker", setdiff(colnames(tbl), "marker"))]
  tbl
}))

write.table(JOURNEY_SAMPLE,
            file      = PATH_JOURNEY_SAMPLE,
            sep       = ";",
            row.names = FALSE,
            quote     = TRUE)

cat(sprintf("  %d samples x %d markers written to:\n  %s\n\n",
            nrow(STAGE_MR[[MARKERS[1]]]$raw),
            length(MARKERS),
            PATH_JOURNEY_SAMPLE))

# --- 4.2  Rarefaction summary -------------------------------------------------

cat("--- Rarefaction summary ------------------------------\n")
cat(sprintf("  %-6s | %10s | %8s | %8s | %12s\n",
            "Marker", "Threshold", "Samples", "ASVs", "Reads"))
cat(sprintf("  %s\n", strrep("-", 54)))

for (marker in MARKERS) {
  mr <- DS[[marker]]$MR_rare
  cat(sprintf("  %-6s | %10s | %8d | %8d | %12s\n",
              marker,
              format(RARE_THRESHOLDS[[marker]], big.mark = ","),
              nrow(mr), ncol(mr),
              format(sum(mr), big.mark = ",")))
}

cat(sprintf("\n  Iterations : %d | Seed : %d\n\n", RARE_ITER, RARE_SEED))

# =============================================================================
# SECTION 5 — RESIDUAL CORRECTION
# =============================================================================
#
# Sequencing depth effect on Observed ASV richness is corrected by extracting
# residuals from a linear model: Observed ~ f(reads), where f is the best-
# performing transformation (linear, log, or sqrt), selected by lowest mean
# |rho| between residuals and read depth across all markers.
#
# Output per marker:
#   $MR_cor  — named numeric vector of residuals (one value per sample)
#   $PARAMS  — all filtering and analytical decisions consolidated
#
# Only indices listed in INDICES_RESIDUAL receive residual correction.
# All other indices use rarefied values from MR_rare directly.
# Reference: Schiaffino et al. (2016) Microbial Ecology 71:819-832
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 5: Residual correction\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Index corrected    :", paste(INDICES_RESIDUAL, collapse = ", "), "\n")
cat("  Significance level :", ALPHA_SIG, "\n")
cat("=============================================================================\n\n")

TRANSFORM_EVAL <- list()
BEST_TRANSFORM <- list()

for (marker in MARKERS) {
  
  cat(sprintf("--- %s -------------------------------------------\n", marker))
  
  mr_clean <- DS[[marker]]$MR_clean
  reads    <- rowSums(mr_clean)
  observed <- vegan::specnumber(round(mr_clean))
  
  # --- 5.1  Evaluate transformations ------------------------------------------
  eval_tbl <- evaluate_transforms(
    values = observed,
    reads  = reads,
    marker = marker,
    index  = "Observed",
    alpha  = ALPHA_SIG
  )
  TRANSFORM_EVAL[[marker]] <- eval_tbl
  best                     <- select_best_transform(eval_tbl)
  BEST_TRANSFORM[[marker]] <- best
  
  # --- 5.2  Console validation table ------------------------------------------
  cat(sprintf("  Best transformation : %s\n\n", best))
  cat(sprintf("  %-8s | %-10s | %-10s | %-10s | %-10s | %-8s | %s\n",
              "Transform", "rho_before", "p_before",
              "rho_after", "p_after", "Corrected", "Validated"))
  cat(sprintf("  %s\n", strrep("-", 75)))
  
  for (i in seq_len(nrow(eval_tbl))) {
    r <- eval_tbl[i, ]
    cat(sprintf("  %-8s | %-10s | %-10s | %-10s | %-10s | %-8s | %s\n",
                r$transform,
                round(r$rho_before, 3), round(r$p_before, 4),
                round(r$rho_after,  3), round(r$p_after,  4),
                ifelse(isTRUE(r$correction_applied), "YES", "NO"),
                ifelse(isTRUE(r$validated),          "YES", "NO [flag]")))
  }
  cat("\n")
  
  # --- 5.3  Apply residual correction -----------------------------------------
  best_row  <- eval_tbl[eval_tbl$transform == best, ]
  mr_cor    <- if (isTRUE(best_row$correction_applied)) {
    cat(sprintf(
      "  Depth effect significant (rho = %s, p = %s) — correction applied.\n",
      round(best_row$rho_before, 3), round(best_row$p_before, 4)))
    setNames(
      compute_residuals(observed, reads, best),
      rownames(mr_clean)
    )
  } else {
    cat(sprintf(
      "  No significant depth effect (rho = %s, p = %s) — raw Observed used.\n",
      round(best_row$rho_before, 3), round(best_row$p_before, 4)))
    setNames(as.numeric(observed), rownames(mr_clean))
  }
  
  DS[[marker]]$MR_cor <- mr_cor
  
  cat(sprintf("  MR_cor : %d values | range [%.2f, %.2f]\n\n",
              length(mr_cor), min(mr_cor, na.rm = TRUE),
              max(mr_cor, na.rm = TRUE)))
  
  # --- 5.4  Build PARAMS ------------------------------------------------------
  DS[[marker]]$PARAMS <- list(
    # Filtering
    exclude_taxa        = EXCLUDE_TAXA[[marker]],
    exclude_samples     = EXCLUDE_SAMPLES[[marker]],
    filter_min_reads    = FILTER_MIN_READS,
    filter_min_samples  = FILTER_MIN_SAMPLES,
    filter_min_identity = FILTER_MIN_IDENTITY,
    filter_min_boot_dom = FILTER_MIN_BOOT_DOMAIN,
    # Rarefaction
    rare_threshold      = RARE_THRESHOLDS[[marker]],
    rare_iter           = RARE_ITER,
    rare_seed           = RARE_SEED,
    # Residual correction
    best_transform      = best,
    indices_residual    = INDICES_RESIDUAL,
    correction_applied  = isTRUE(best_row$correction_applied),
    correction_validated = isTRUE(best_row$validated),
    transform_eval      = eval_tbl,
    # Journey statistics
    n_samples_raw       = nrow(DS[[marker]]$MR),
    n_asvs_raw          = ncol(DS[[marker]]$MR),
    n_samples_clean     = nrow(DS[[marker]]$MR_clean),
    n_asvs_clean        = ncol(DS[[marker]]$MR_clean),
    n_samples_rare      = nrow(DS[[marker]]$MR_rare),
    n_asvs_rare         = ncol(DS[[marker]]$MR_rare)
  )
  
}  # end marker loop

# --- 5.5  Cross-marker summary ------------------------------------------------

cat("--- Transformation summary ---------------------------\n")
cat(sprintf("  %-6s | %-10s | %-10s | %-10s | %-8s | %s\n",
            "Marker", "Transform", "rho_before", "rho_after",
            "Corrected", "Validated"))
cat(sprintf("  %s\n", strrep("-", 66)))

for (marker in MARKERS) {
  best <- BEST_TRANSFORM[[marker]]
  r    <- TRANSFORM_EVAL[[marker]][
    TRANSFORM_EVAL[[marker]]$transform == best, ]
  cat(sprintf("  %-6s | %-10s | %-10s | %-10s | %-8s | %s\n",
              marker, best,
              round(r$rho_before, 3), round(r$rho_after, 3),
              ifelse(isTRUE(r$correction_applied), "YES", "NO"),
              ifelse(isTRUE(r$validated),          "YES", "NO [flag]")))
}
cat("\n")

# =============================================================================
# SECTION 5.5 — FINAL VALIDATION
# =============================================================================
#
# Pre-save validation confirming that filtering, rarefaction and residual
# correction performed as expected. Two diagnostics:
#
#   DS_01_55_reads_vs_asvs : reads vs ASVs at three pipeline stages
#                               (MR_clean, MR_rare, MR_cor)
#   DS_01_55_design_balance: sample counts per treatment x timepoint x block
#                               at MR_clean and MR_rare levels
#
# x-axis for all scatter panels = raw reads (DS[[marker]]$MR) matched
# by rowname to samples present in each stage matrix.
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 5.5: Final validation\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# =============================================================================
# 5.5.1 — READS vs ASVs AT THREE PIPELINE STAGES
# =============================================================================
# Row per marker, three columns:
#   Col 1 — MR_clean  : x = raw reads, y = Observed ASVs from MR_clean
#   Col 2 — MR_rare   : x = raw reads, y = Observed ASVs from MR_rare
#   Col 3 — MR_cor    : x = raw reads, y = log-residuals from MR_cor
# Spearman rho per panel. Red = significant (p <= ALPHA_SIG).
# =============================================================================

cat("--- 5.5.1  Reads vs ASVs across pipeline stages -----\n\n")

#' Build one scatter panel: raw reads (x) vs diversity metric (y).
#'
#' x is always derived from the raw MR matrix, matched by rowname to the
#' samples present in the stage matrix. This ensures honest reporting of
#' the original sequencing depth even when samples have been excluded.
#'
#' @param raw_reads  Named numeric. rowSums of raw MR (all samples).
#' @param y_values   Named numeric. Metric values (samples in stage matrix).
#' @param color      Character. Point and line colour.
#' @param title      Character. Panel title.
#' @param y_lab      Character. y-axis label.
#' @param alpha_sig  Numeric. Significance threshold for Spearman p-value.
#' @return A ggplot object.
plot_reads_vs_metric <- function(raw_reads, y_values,
                                 color, title, y_lab,
                                 alpha_sig = 0.05) {
  
  common <- intersect(names(raw_reads), names(y_values))
  df     <- data.frame(
    reads  = raw_reads[common],
    metric = y_values[common]
  )
  
  test    <- cor.test(df$reads, df$metric, method = "spearman", exact = FALSE)
  rho     <- round(test$estimate, 3)
  p_val   <- test$p.value
  p_label <- ifelse(p_val < 0.001, "p < 0.001",
                    paste0("p = ", round(p_val, 3)))
  sig     <- p_val <= alpha_sig
  
  ggplot(df, aes(x = reads, y = metric)) +
    geom_point(color = color, alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE,
                color     = ifelse(sig, "#CC0000", "grey40"),
                linewidth = 0.8) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",")) +
    labs(title    = title,
         subtitle = paste0("Spearman rho = ", rho, " | ", p_label,
                           " | n = ", nrow(df)),
         x        = "Raw reads per sample",
         y        = y_lab) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(
            size  = 8,
            color = ifelse(sig, "#CC0000", "grey40")))
}

fig_stages_list <- setNames(lapply(MARKERS, function(marker) {
  
  raw_reads <- rowSums(DS[[marker]]$MR)
  
  # Col 1: MR_clean — Observed ASVs
  asvs_clean <- setNames(
    vegan::specnumber(round(DS[[marker]]$MR_clean)),
    rownames(DS[[marker]]$MR_clean)
  )
  p1 <- plot_reads_vs_metric(
    raw_reads = raw_reads,
    y_values  = asvs_clean,
    color     = MARKER_COLORS[marker],
    title     = paste0(MARKER_LABELS[marker], " — MR_clean"),
    y_lab     = "Observed ASVs",
    alpha_sig = ALPHA_SIG
  )
  
  # Col 2: MR_rare — Observed ASVs
  asvs_rare <- setNames(
    vegan::specnumber(round(DS[[marker]]$MR_rare)),
    rownames(DS[[marker]]$MR_rare)
  )
  p2 <- plot_reads_vs_metric(
    raw_reads = raw_reads,
    y_values  = asvs_rare,
    color     = MARKER_COLORS[marker],
    title     = paste0(MARKER_LABELS[marker], " — MR_rare"),
    y_lab     = "Observed ASVs",
    alpha_sig = ALPHA_SIG
  )
  
  # Col 3: MR_cor — log-residuals
  p3 <- plot_reads_vs_metric(
    raw_reads = raw_reads,
    y_values  = DS[[marker]]$MR_cor,
    color     = MARKER_COLORS[marker],
    title     = paste0(MARKER_LABELS[marker], " — MR_cor"),
    y_lab     = "Observed ASVs (log-residuals)",
    alpha_sig = ALPHA_SIG
  )
  
  # Console summary
  cat(sprintf("  %s:\n", marker))
  for (stage in c("MR_clean", "MR_rare", "MR_cor")) {
    y <- if (stage == "MR_clean") asvs_clean else
      if (stage == "MR_rare")  asvs_rare  else
        DS[[marker]]$MR_cor
    common <- intersect(names(raw_reads), names(y))
    test   <- cor.test(raw_reads[common], y[common],
                       method = "spearman", exact = FALSE)
    cat(sprintf("    %-10s : rho = %6.3f | p = %.4f | n = %d | %s\n",
                stage,
                round(test$estimate, 3),
                test$p.value,
                length(common),
                ifelse(test$p.value <= ALPHA_SIG,
                       "significant [flag]", "OK")))
  }
  cat("\n")
  
  wrap_plots(p1, p2, p3, nrow = 1)
  
}), MARKERS)

fig_stages_combined <- wrap_plots(fig_stages_list, ncol = 1) +
  plot_annotation(
    title    = paste0(PROJECT_NAME,
                      " — Reads vs ASVs across pipeline stages"),
    subtitle = paste0(
      "x = raw reads (before any filter) | ",
      "Col 1: MR_clean | Col 2: MR_rare | Col 3: MR_cor (log-residuals) | ",
      "Red = significant depth effect (p <= ", ALPHA_SIG, ")"),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

print(fig_stages_combined)
save_figure(fig_stages_combined, "01_55_reads_vs_asvs",
            width = 14, height = 5 * length(MARKERS))


cat("--- Section 5.5 complete — proceed to Section 6 -----\n\n")

# =============================================================================
# SECTION 6 — SAVE DS.RData
# =============================================================================

cat("--- Saving", paste0(PROJECT_NAME, ".RData"), "-------------------\n\n")

save(DS, file = PATH_RDATA)

cat(sprintf("  Saved: %s\n\n", PATH_RDATA))

cat("  Object structure:\n")
cat(sprintf("  %-6s | %-12s | %-8s | %-8s | %s\n",
            "Marker", "Slot", "Samples", "ASVs", "Reads"))
cat(sprintf("  %s\n", strrep("-", 54)))

for (marker in MARKERS) {
  for (slot in c("MR", "MR_clean", "MR_rare")) {
    mr <- DS[[marker]][[slot]]
    cat(sprintf("  %-6s | %-12s | %8d | %8d | %12s\n",
                marker, slot,
                nrow(mr), ncol(mr),
                format(sum(mr), big.mark = ",")))
  }
  cat(sprintf("  %-6s | %-12s | %8d values | range [%.2f, %.2f]\n",
              marker, "MR_cor",
              length(DS[[marker]]$MR_cor),
              min(DS[[marker]]$MR_cor, na.rm = TRUE),
              max(DS[[marker]]$MR_cor, na.rm = TRUE)))
  cat(sprintf("  %s\n", strrep("-", 54)))
}

cat("\n  Next step: source('02_dataset_quality.R')\n\n")

# =============================================================================
# SECTION 7 — SESSION INFO
# =============================================================================

cat("--- Saving session info ------------------------------\n")
sink(PATH_SESSION_TXT)
cat(PROJECT_NAME, "— 01_filtering.R session information\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
cat("  Saved:", PATH_SESSION_TXT, "\n\n")


# =============================================================================
# SECTION 8 — CLEAN ENVIRONMENT
# =============================================================================

rm(list = setdiff(ls(), c(
  "DS",
  "ROOT_DIR",
  "PROJECT_NAME",
  "MARKERS"
)))

cat("  Environment cleaned.\n")
cat("  Objects retained: DS, ROOT_DIR, PROJECT_NAME, MARKERS\n")
cat("\n  Next step: source('02_dataset_quality.R')\n")
cat("=============================================================================\n")

