# =============================================================================
# File      : 00_import.R
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Author    : David Singer
# License   : MIT
# =============================================================================
#
# DESCRIPTION
# -----------
# First script of the MetaBarFlow statistics stage. Imports the pipeline's
# per-marker outputs (read matrix + taxonomic assignment table) and saves them
# as a single standardised R object. No matrix modification is performed here.
#
# PIPELINE CONTEXT
# ----------------
# Raw sequencing data are processed upstream by the MetaBarFlow bioinformatics
# pipeline (scripts/): merge -> quality check -> primer trimming (Cutadapt) ->
# ASV inference + chimera removal (DADA2) -> taxonomic assignment (VSEARCH +
# DADA2 assignTaxonomy). This script consumes the outputs of steps 4 and 5.
#
# INPUT FILES (read from the pipeline's results/ directory)
# ---------------------------------------------------------
# For each marker (e.g. 18S, ITS, 16S), two files are required:
#
# 1. MR file (read matrix) — results/<MARKER>/04_ASV/<MARKER>_MR.csv
#    Rows    : samples (sample ID as row name, first column, no header)
#    Columns : ASVs (ASV_ID as column header)
#    Values  : raw read counts (integers >= 0)
#    Example :
#      ;ASV_001;ASV_002;...
#      S1;152;0;...
#      S2;0;93;...
#
# 2. ASS file (taxonomic assignment) — results/<MARKER>/05_taxonomy/<MARKER>_ASS.csv
#    Rows    : ASVs (one per row)
#    Columns : ASV_ID, Nb_reads, Pct_identity, GenBank, Domain, Supergroup,
#              Division, Subdivision, Class, Order, Family, Genus, Species,
#              <rank>_boot (9 bootstrap columns, values 0–1), Sequence
#    Constraint: ASV set must be identical to MR columns, in the same order.
#
# OUTPUT FILES
# ------------
# <PROJECT>.RData                      : raw import object (all markers)
# outputs/<PROJECT>_sessioninfo_00.txt : R session information
#
# OBJECT STRUCTURE
# ----------------
# DS
# ├── [["18S"]]
# │   ├── $MR   — raw read matrix (samples x ASVs)
# │   └── $ASS  — full taxonomic assignment table (ASVs x fields)
# ├── [["ITS"]] — same structure
# └── [["16S"]] — same structure
# Slots $MR_clean, $ASS_clean, $MR_rare, $ASS_rare, $MR_cor, $PARAMS are added
# by 01_filtering.R.
#
# USAGE
# -----
# 1. Run the bioinformatics pipeline first so results/<MARKER>/ exists.
# 2. Open Stats.Rproj in RStudio (sets the working directory to Stats/),
#    or setwd() to the Stats/ folder. Verify with getwd().
# 3. Edit Section 1 parameters if needed (markers, paths, delimiters).
# 4. source("00_import.R")  — or: Rscript 00_import.R
# 5. Continue with 01_filtering.R
#
# DEPENDENCIES
# ------------
# R >= 4.0.0. No external packages required for this script.
# =============================================================================


# =============================================================================
# SECTION 1 — PARAMETERS
# =============================================================================

# --- 1.1  Project name --------------------------------------------------------
# Names the output object (.RData) and all output files. Change per project.

PROJECT_NAME <- "MetaBarFlow"

# --- 1.2  Paths ---------------------------------------------------------------
# Paths are built relative to the Stats/ working directory. Open Stats.Rproj in
# RStudio (working directory set automatically) or setwd() to Stats/.
# RESULTS_DIR points at the pipeline output (../results by default). Edit it if
# your results live elsewhere.

ROOT_DIR    <- getwd()
RESULTS_DIR <- normalizePath(file.path(ROOT_DIR, "..", "results"), mustWork = FALSE)

# Create standard output folder structure if not already present
for (d in file.path(ROOT_DIR, "outputs", c("figures", "tables")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

PATH_RDATA       <- file.path(ROOT_DIR, paste0(PROJECT_NAME, ".RData"))
PATH_SESSION_TXT <- file.path(ROOT_DIR, "outputs",
                              paste0(PROJECT_NAME, "_sessioninfo_00.txt"))

# --- 1.3  Markers to process --------------------------------------------------
# List the markers to import. Remove any marker not available.

MARKERS <- c("18S", "ITS", "16S")

# --- 1.4  Marker labels and colours -------------------------------------------
# Used consistently across all figures. Edit to match your markers.

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

# --- 1.5  Input file paths ----------------------------------------------------
# Built automatically from MARKERS and RESULTS_DIR following the pipeline's
# output layout. No manual editing needed unless your layout differs.

PATH_MR <- setNames(
  lapply(MARKERS, function(m)
    file.path(RESULTS_DIR, m, "04_ASV", paste0(m, "_MR.csv"))),
  MARKERS)
PATH_ASS <- setNames(
  lapply(MARKERS, function(m)
    file.path(RESULTS_DIR, m, "05_taxonomy", paste0(m, "_ASS.csv"))),
  MARKERS)

# --- 1.6  File delimiters -----------------------------------------------------
# The pipeline writes both files semicolon-delimited. Adjust if you changed it.

DELIMITER_MR  <- ";"   # MR file separator
DELIMITER_ASS <- ";"   # ASS file separator

# --- 1.7  ASS column names ----------------------------------------------------
# Names of key columns in the ASS file (as written by 05_taxonomy.R).
# Edit only if you renamed columns.

ASS_COL_ID       <- "ASV_ID"       # ASV identifier column
ASS_COL_IDENTITY <- "Pct_identity" # VSEARCH % identity column
ASS_COL_BOOT_DOM <- "Domain_boot"  # bootstrap support at Domain level
ASS_COL_DIVISION <- "Division"     # taxonomic Division column
ASS_COL_SUBDIV   <- "Subdivision"  # taxonomic Subdivision column

# =============================================================================
# SECTION 2 — IMPORT FUNCTIONS
# =============================================================================

# --- 2.1  Pre-flight file existence check -------------------------------------

#' Verify that all required input files exist before starting import.
#'
#' Checks all MR and ASS file paths. Stops execution with an informative error
#' listing all missing files if any are not found.
#'
#' @param path_mr   Named list. MR file paths keyed by marker name.
#' @param path_ass  Named list. ASS file paths keyed by marker name.
#' @param markers   Character vector. Markers to check.
preflight_check <- function(path_mr, path_ass, markers) {
  all_paths <- c(
    setNames(unlist(path_mr[markers]),  paste(markers, "MR")),
    setNames(unlist(path_ass[markers]), paste(markers, "ASS"))
  )
  missing <- all_paths[!file.exists(all_paths)]
  if (length(missing) > 0)
    stop("\nPre-flight check failed — missing input files:\n",
         paste0("  [", names(missing), "] ", missing, collapse = "\n"),
         "\nHas the bioinformatics pipeline been run? Check RESULTS_DIR ",
         "and MARKERS in Section 1.")
  cat("  Pre-flight check  : all input files found [OK]\n")
}

# --- 2.2  ASV matrix import ---------------------------------------------------

#' Import a community matrix (MR) from a delimited CSV file.
#'
#' Reads a semicolon-delimited (or configurable) file with samples as rows
#' and ASVs as columns. The first column is used as row names (sample IDs).
#' Values are expected to be raw read counts (integers >= 0).
#'
#' @param path      Character. Path to the MR CSV file.
#' @param marker    Character. Marker name for console output.
#' @param delimiter Character. Field separator.
#' @return A numeric data.frame (samples x ASVs) with sample IDs as rownames.
import_mr <- function(path, marker, delimiter = ";") {
  cat("  Importing", marker, "MR ...\n")
  mr <- read.table(path,
                   sep         = delimiter,
                   header      = TRUE,
                   row.names   = 1,
                   check.names = FALSE)
  # Validate read counts — must be non-negative integers
  if (any(is.na(mr)))
    stop(marker, " MR: NA values detected — check delimiter and file format")
  if (any(mr < 0, na.rm = TRUE))
    stop(marker, " MR: negative values detected — raw read counts must be >= 0")
  if (any(mr != round(mr), na.rm = TRUE))
    warning(marker, " MR: non-integer values detected — verify input file")
  if (any(duplicated(rownames(mr))))
    stop(marker, " MR: duplicate sample IDs detected — check input file")

  cat("    Samples :", nrow(mr), "\n")
  cat("    ASVs    :", ncol(mr), "\n")
  cat("    Reads   :", format(sum(mr), big.mark = ","), "\n")
  return(mr)
}

# --- 2.3  Taxonomic assignment table import -----------------------------------

#' Import a taxonomic assignment table (ASS) from a delimited CSV file.
#'
#' Reads a semicolon-delimited (or configurable) file with one ASV per row.
#' The column identified by col_id is used as row names (ASV IDs).
#'
#' @param path      Character. Path to the ASS CSV file.
#' @param marker    Character. Marker name for console output.
#' @param col_id    Character. Name of the ASV identifier column.
#' @param delimiter Character. Field separator.
#' @return A data.frame (ASVs x fields) with ASV IDs as rownames.
import_ass <- function(path, marker, col_id = "ASV_ID", delimiter = ";") {
  cat("  Importing", marker, "ASS ...\n")
  ass <- read.table(path,
                    sep              = delimiter,
                    header           = TRUE,
                    row.names        = NULL,
                    check.names      = FALSE,
                    stringsAsFactors = FALSE,
                    quote            = "\"")
  if (!col_id %in% colnames(ass))
    stop(marker, " ASS: column '", col_id, "' not found. ",
         "Check ASS_COL_ID in Section 1.7.")
  if (any(duplicated(ass[[col_id]])))
    stop(marker, " ASS: duplicate ASV IDs in column '", col_id,
         "' — check input file")
  rownames(ass) <- ass[[col_id]]
  cat("    ASVs    :", nrow(ass), "\n")
  return(ass)
}

# --- 2.4  Consistency check ---------------------------------------------------

#' Check that MR column names and ASS row names are identical (same set,
#' same order). Stops with an informative error if a mismatch is detected.
#'
#' Called immediately after raw import, before any matrix manipulation.
#' Identity is required: every ASV in the MR must have an assignment row,
#' and no extra rows are permitted in the ASS.
#'
#' @param mr     data.frame. Raw samples x ASVs matrix.
#' @param ass    data.frame. Raw taxonomic assignment table.
#' @param marker Character. Marker name for the error message.
check_mr_ass <- function(mr, ass, marker) {
  if (identical(colnames(mr), rownames(ass))) {
    cat("    MR / ASS consistency : OK (", ncol(mr), "ASVs matched)\n")
  } else {
    n_mr  <- ncol(mr)
    n_ass <- nrow(ass)
    only_mr  <- setdiff(colnames(mr), rownames(ass))
    only_ass <- setdiff(rownames(ass), colnames(mr))
    stop(marker, ": MR / ASS mismatch.\n",
         "  MR columns : ", n_mr,  " ASVs\n",
         "  ASS rows   : ", n_ass, " ASVs\n",
         if (length(only_mr)  > 0) paste0("  In MR not ASS: ",
                                          paste(head(only_mr,  5), collapse = "; "),
                                          if (length(only_mr)  > 5) " ...", "\n") else "",
         if (length(only_ass) > 0) paste0("  In ASS not MR: ",
                                          paste(head(only_ass, 5), collapse = "; "),
                                          if (length(only_ass) > 5) " ...", "\n") else "",
         "Check pipeline outputs for ", marker, ".")
  }
}

# =============================================================================
# SECTION 3 — IMPORT
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Data import\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Stats root   :", ROOT_DIR, "\n")
cat("  Results dir  :", RESULTS_DIR, "\n")
cat("  Markers      :", paste(MARKERS, collapse = "; "), "\n")
cat("=============================================================================\n\n")

# --- 3.1  Pre-flight check ----------------------------------------------------

cat("--- Pre-flight check ----------------------------------\n")
preflight_check(PATH_MR, PATH_ASS, MARKERS)

# --- 3.2  MR and ASS per marker -----------------------------------------------

DS <- list()

for (marker in MARKERS) {

  cat("\n---", marker, "----------------------------------------------\n")

  mr  <- import_mr(PATH_MR[[marker]],  marker, delimiter = DELIMITER_MR)
  ass <- import_ass(PATH_ASS[[marker]], marker,
                    col_id    = ASS_COL_ID,
                    delimiter = DELIMITER_ASS)

  DS[[marker]] <- list(MR = mr, ASS = ass)

}  # end marker loop

# =============================================================================
# SECTION 4 — CONSISTENCY CHECKS
# =============================================================================

cat("\n--- Consistency checks --------------------------------\n")

# Note: check_mr_ass() calls stop() on failure. If a check fails mid-loop, DS is
# partially built in memory but <PROJECT>.RData is never written (save() occurs
# in Section 6 only). Re-run after fixing the input files — no partial save risk.
for (marker in MARKERS) {
  cat("\n  [", marker, "]\n", sep = "")
  check_mr_ass(DS[[marker]]$MR, DS[[marker]]$ASS, marker)
}

cat("\n  All consistency checks passed.\n")

# --- Zero-read sample check ---------------------------------------------------
# Samples with zero reads across all ASVs will cause division-by-zero errors
# in relative abundance calculations downstream. Flagged here as warnings —
# consider adding to EXCLUDE_SAMPLES in 01_filtering.R.

cat("\n--- Zero-read sample check ----------------------------\n")
for (marker in MARKERS) {
  zero_samp <- rownames(DS[[marker]]$MR)[rowSums(DS[[marker]]$MR) == 0]
  if (length(zero_samp) > 0)
    warning(marker, ": ", length(zero_samp),
            " sample(s) with zero reads — ",
            paste(zero_samp, collapse = "; "))
  else
    cat(sprintf("  %-6s | no zero-read samples [OK]\n", marker))
}

# =============================================================================
# SECTION 5 — SUMMARY PER MARKER
# =============================================================================

cat("\n--- Raw data summary ----------------------------------\n")
cat(sprintf("  %-6s | %7s | %8s | %14s | %10s | %10s | %10s\n",
            "Marker", "Samples", "ASVs", "Reads",
            "Mean/samp", "Min/samp", "Max/samp"))
cat(sprintf("  %s\n", strrep("-", 72)))

for (marker in MARKERS) {
  mr <- DS[[marker]]$MR
  rs <- rowSums(mr)
  cat(sprintf("  %-6s | %7d | %8d | %14s | %10s | %10s | %10s\n",
              marker,
              nrow(mr),
              ncol(mr),
              format(sum(mr),         big.mark = ","),
              format(round(mean(rs)), big.mark = ","),
              format(min(rs),         big.mark = ","),
              format(max(rs),         big.mark = ",")))
}
cat(sprintf("  %s\n", strrep("-", 72)))
cat("  Note: values are raw counts — no filtering applied.\n")
cat("  Note: any taxonomic filter (e.g. non-target removal for 18S) is applied\n")
cat("        in 01_filtering.R via EXCLUDE_TAXA.\n")

# =============================================================================
# SECTION 6 — SAVE <PROJECT>.RData
# =============================================================================

cat("\n--- Saving", paste0(PROJECT_NAME, ".RData"), "-------------------\n")
save(DS, file = PATH_RDATA)
cat("  Saved:", PATH_RDATA, "\n")
cat("  Object: DS —", length(MARKERS), "markers\n")
cat("  Slots per marker: $MR, $ASS\n")
cat("  Note: $MR_clean, $ASS_clean, $MR_rare, $MR_cor, $PARAMS added by 01_filtering.R\n")

# =============================================================================
# SECTION 7 — SESSION INFO
# =============================================================================

cat("\n--- Saving session info ------------------------------\n")
sink(PATH_SESSION_TXT)
cat(PROJECT_NAME, "— 00_import.R session information\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
cat("  Saved:", PATH_SESSION_TXT, "\n")

# =============================================================================
# SECTION 8 — CLEAN ENVIRONMENT
# =============================================================================

rm(list = setdiff(ls(), c(
  "DS",
  "ROOT_DIR",
  "PROJECT_NAME",
  "MARKERS",
  "MARKER_LABELS",
  "MARKER_COLORS"
)))

cat("\n  Environment cleaned.\n")
cat("  Objects retained: DS, ROOT_DIR, PROJECT_NAME, MARKERS,\n")
cat("                    MARKER_LABELS, MARKER_COLORS\n")
cat("=============================================================================\n")
