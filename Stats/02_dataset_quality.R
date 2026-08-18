# =============================================================================
# File      : 02_dataset_quality.R
# Project   : DS – Drought Research at Yvorne 2050
# Author    : David Singer
# Institution: HES-SO Changins, Nyon, Switzerland
# Date      : 2026-04-28
# Version   : 1.0
# License   : MIT
# =============================================================================
#
# DESCRIPTION
# -----------
# Dataset quality control and diagnostic visualisation for the DS
# statistical pipeline. Loads DS.RData produced by 01_filtering.R
# and generates a comprehensive set of diagnostic figures and tables.
# No matrix modification is performed — all outputs are figures and tables only.
#
# PIPELINE CONTEXT
# ----------------
# This script is the third in the MarkerFlow statistics stage:
#   00_import.R    → imports pipeline outputs → <PROJECT>.RData ($MR, $ASS)
#   01_filtering.R → filters and corrects → <PROJECT>.RData ($MR_clean,
#                    $MR_rare, $MR_cor, $PARAMS)
#   02_dataset_quality.R → diagnostic figures and tables (this script)
#
# INPUT
# -----
# <PROJECT>.RData — produced by 01_filtering.R
# Contains per marker: $MR, $MR_clean, $MR_rare, $MR_cor,
#                      $ASS, $ASS_clean, $ASS_rare, $PARAMS
# This stage is metadata-free.
#
# OUTPUT
# ------
# outputs/figures/02_*.pdf / *.jpg : diagnostic figures
# outputs/tables/02_*.csv          : diagnostic tables
# outputs/[PROJECT]_sessioninfo_02.txt : session information
# No new RData produced.
#
# SECTIONS
# --------
#  1. Parameters
#  LOAD DATA
#  2. Read depth distribution
#  3. Library size vs richness (raw reads vs final ASVs — all pipeline stages)
#  4. Rarefaction curves
#  5. ASV sequence length distribution
#  6. Identity distribution (Mahé et al. 2017)
#  7. Filtering journey visualisation
#  8. Session info
#  9. Clean environment
#
# DEPENDENCIES
# ------------
# ggplot2, patchwork, scales, vegan
#
# REFERENCES
# ----------
# Mahé et al. (2017) Swarm v2. PeerJ 5:e3166
# Oksanen et al. (2024) vegan. R package.
# =============================================================================


# =============================================================================
# SECTION 1 — PARAMETERS
# =============================================================================

# --- 1.1  Project identifiers -------------------------------------------------

PROJECT_NAME <- "MarkerFlow"
MARKERS      <- c("18S", "ITS", "16S")

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

# --- 1.2  Paths ---------------------------------------------------------------
# All paths built relative to the working directory (project root).
# Open the .Rproj file in RStudio — working directory is set automatically.

ROOT_DIR <- getwd()

PATH_RDATA       <- file.path(ROOT_DIR, paste0(PROJECT_NAME, ".RData"))
PATH_FIG_DIR     <- file.path(ROOT_DIR, "outputs", "figures")
PATH_TAB_DIR     <- file.path(ROOT_DIR, "outputs", "tables")
PATH_SESSION_TXT <- file.path(ROOT_DIR, "outputs",
                              paste0(PROJECT_NAME, "_sessioninfo_02.txt"))

dir.create(PATH_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PATH_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1.3  Expected amplicon length ranges (bp) --------------------------------
# Used to flag ASVs outside the expected size range in Section 5.
# Set to NULL for a marker to skip the length filter flag.

AMPLICON_LENGTHS <- list(
  "18S" = c(350, 500),
  "ITS" = c(200, 600),
  "16S" = c(300, 550)
)

# --- 1.5  Rarefaction parameters (from 01_filtering.R) -----------------------
# Used to overlay thresholds on rarefaction curves (Section 4).
# These must match the values used in 01_filtering.R.
# Retrieved automatically from DS$PARAMS after loading — see LOAD DATA.

# --- 1.6  Identity distribution parameters ------------------------------------
# Bin width for % identity histograms (Section 6).

IDENTITY_BINWIDTH <- 0.5   # percentage points

# --- 1.7  Statistical parameters ----------------------------------------------
# Significance threshold for Spearman correlation annotations (Section 3).

ALPHA_SIG <- 0.05

# --- 1.8  Rarefaction curve step sizes ----------------------------------------
# Step size for rarefaction curve interpolation per marker (Section 4).
# Smaller = smoother curves but slower. Increase for large datasets.

RARE_STEP <- list(
  "18S" = 100,
  "ITS" = 200,
  "16S" = 500
)

# --- 1.9  Figure options ------------------------------------------------------

SAVE_PDF        <- TRUE    # save figures as PDF
SAVE_JPEG       <- TRUE    # save figures as JPEG
FIG_DPI         <- 300     # resolution for JPEG output
COMBINE_MARKERS <- TRUE    # TRUE  = one multi-panel figure per analysis
                           # FALSE = one figure per marker

# =============================================================================
# LOAD DATA
# =============================================================================
# Loads DS.RData produced by 01_filtering.R.
# Retrieves rarefaction thresholds from DS$PARAMS for use in figures.
# No matrix modification performed here.
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— 02_dataset_quality.R\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Loading", PATH_RDATA, "...\n")
cat("=============================================================================\n\n")

if (!file.exists(PATH_RDATA))
  stop("Input file not found: ", PATH_RDATA,
       "\nRun 01_filtering.R first or check ROOT_DIR.")

load(PATH_RDATA)
cat("  DS loaded —", length(MARKERS), "markers:",
    paste(MARKERS, collapse = ", "), "\n")

# --- Retrieve rarefaction thresholds from PARAMS ------------------------------
# PARAMS are stored per marker by 01_filtering.R.
# Retrieved here for use in rarefaction curve overlays (Section 4).

RARE_THRESHOLDS <- setNames(
  lapply(MARKERS, function(m) {
    t <- DS[[m]]$PARAMS$rare_threshold
    if (is.null(t)) {
      warning(m, ": PARAMS$rare_threshold not found — ",
              "rarefaction threshold overlay will be skipped.")
      return(NA)
    }
    return(t)
  }),
  MARKERS
)

cat("\n  Rarefaction thresholds (from DS$PARAMS):\n")
for (m in MARKERS)
  cat(sprintf("    %-6s : %s reads\n",
              m, ifelse(is.na(RARE_THRESHOLDS[[m]]),
                        "not found",
                        format(RARE_THRESHOLDS[[m]], big.mark = ","))))

# --- Quick slot availability check -------------------------------------------
# Verify expected slots are present for each marker.
# Downstream sections will fail gracefully if a slot is missing.

cat("\n  Slot availability per marker:\n")
expected_slots <- c("MR", "MR_clean", "MR_rare", "MR_cor",
                    "ASS", "ASS_clean", "ASS_rare")
for (m in MARKERS) {
  present <- expected_slots[expected_slots %in% names(DS[[m]])]
  missing <- setdiff(expected_slots, names(DS[[m]]))
  cat(sprintf("    %-6s | present: %-40s",
              m, paste(present, collapse = ", ")))
  if (length(missing) > 0)
    cat(sprintf(" | MISSING: %s", paste(missing, collapse = ", ")))
  cat("\n")
}
cat("\n")

# --- Shared utility functions -------------------------------------------------
# Defined here once — available to all sections.

#' Save a ggplot figure as PDF and/or JPEG depending on SAVE_PDF/SAVE_JPEG flags.
#'
#' @param fig    ggplot object. Figure to save.
#' @param name   Character. Output filename without extension.
#' @param width  Numeric. Figure width in inches.
#' @param height Numeric. Figure height in inches.
save_figure <- function(fig, name, width, height) {
  if (SAVE_PDF)
    ggsave(file.path(PATH_FIG_DIR, paste0(name, ".pdf")),
           fig, width = width, height = height)
  if (SAVE_JPEG)
    ggsave(file.path(PATH_FIG_DIR, paste0(name, ".jpg")),
           fig, width = width, height = height, dpi = FIG_DPI)
}

# =============================================================================
# SECTION 2 — READ DEPTH DISTRIBUTION
# =============================================================================
# Visualises sequencing depth across samples at the raw stage and shows
# the impact of rarefaction on sample retention.
#
# Figures produced:
#   2.1 — Violin + jitter plot of raw read depth per marker
#          Subtitle: n | median | min | max
#   2.2 — Bar chart of samples retained vs excluded by rarefaction
#          Based on MR_clean (filtered matrix) vs rarefaction threshold
#
# Summary table saved as CSV.
#
# Scientific rationale:
# Figure 2.1 uses raw reads ($MR) as the basis — this is the true sequencing
# effort before any processing and provides the most honest representation
# of depth variation across samples. Filtered reads ($MR_clean) would
# underestimate the original depth variation.
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 2: Read depth distribution\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
library(ggplot2)
library(patchwork)
library(scales)

# --- Build read depth summary per marker --------------------------------------

depth_stats <- setNames(lapply(MARKERS, function(m) {
  rs <- rowSums(DS[[m]]$MR)
  list(
    n      = length(rs),
    median = median(rs),
    mean   = round(mean(rs)),
    min    = min(rs),
    max    = max(rs),
    reads  = rs
  )
}), MARKERS)

# Console summary
cat("--- Raw read depth summary ---------------------------\n\n")
cat(sprintf("  %-6s | %5s | %10s | %10s | %10s\n",
            "Marker", "n", "Median", "Min", "Max"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (m in MARKERS) {
  s <- depth_stats[[m]]
  cat(sprintf("  %-6s | %5d | %10s | %10s | %10s\n",
              m, s$n,
              format(s$median, big.mark = ","),
              format(s$min,    big.mark = ","),
              format(s$max,    big.mark = ",")))
}
cat("\n")

# Save summary table
depth_summary <- do.call(rbind, lapply(MARKERS, function(m) {
  s <- depth_stats[[m]]
  data.frame(marker = m, n = s$n, median = s$median,
             mean = s$mean, min = s$min, max = s$max,
             stringsAsFactors = FALSE)
}))
write.table(depth_summary,
            file.path(PATH_TAB_DIR,
                      paste0(PROJECT_NAME, "_02_read_depth_summary.csv")),
            sep = ";", row.names = FALSE, quote = TRUE)
cat("  Summary table saved to:", PATH_TAB_DIR, "\n\n")

# =============================================================================
# FIGURE 2.1 — RAW READ DEPTH VIOLIN + JITTER PER MARKER
# =============================================================================

cat("--- Figure 2.1: Raw read depth violin ---------------\n")

fig2_1_list <- lapply(MARKERS, function(m) {
  
  s   <- depth_stats[[m]]
  col <- MARKER_COLORS[m]
  
  df <- data.frame(
    marker = MARKER_LABELS[m],
    reads  = s$reads,
    stringsAsFactors = FALSE
  )
  
  subtitle <- sprintf("n=%d  |  median=%s  |  min=%s  |  max=%s",
                      s$n,
                      format(s$median, big.mark = ","),
                      format(s$min,    big.mark = ","),
                      format(s$max,    big.mark = ","))
  
  ggplot(df, aes(x = marker, y = reads)) +
    geom_violin(fill = col, color = col,
                alpha = 0.3, linewidth = 0.7) +
    geom_jitter(color = col, alpha = 0.5,
                width = 0.15, size = 1.5) +
    stat_summary(fun = median, geom = "crossbar",
                 width = 0.4, linewidth = 0.7,
                 color = "grey20") +
    scale_y_continuous(labels = scales::comma) +
    labs(title    = MARKER_LABELS[m],
         subtitle = subtitle,
         x        = NULL,
         y        = "Reads per sample") +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(size = 10, face = "bold", color = col),
      plot.subtitle    = element_text(size = 8,  color = "grey40"),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
})

fig2_1 <- wrap_plots(fig2_1_list, nrow = 1) +
  plot_annotation(
    title    = paste0(PROJECT_NAME, " — Raw read depth distribution"),
    subtitle = "Horizontal bar = median | Raw reads before any filtering",
    theme    = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

print(fig2_1)
save_figure(fig2_1, "02_01_read_depth_violin",
            width = 10, height = 6)
cat("  Figure 2.1 saved.\n\n")

# =============================================================================
# SECTION 3 — LIBRARY SIZE VS RICHNESS
# =============================================================================
# Scatter plots showing the relationship between raw sequencing depth and
# Observed ASV richness at each pipeline stage.
#
# x-axis (all panels) : raw reads per sample ($MR) — total before any filter
# y-axis col 1        : Observed ASVs from $MR_clean (after all filters)
# y-axis col 2        : Observed ASVs from $MR_rare (after rarefaction)
# y-axis col 3        : log-residuals from $MR_cor (after correction)
#
# One row per marker, three columns.
# Spearman rho annotated per panel (red = significant at ALPHA_SIG).
# Linear regression line overlaid for visual trend.
#
# Scientific rationale:
# x = raw reads reflects the true sequencing effort before any processing.
# This is the honest measure of depth effect — using filtered or rarefied
# reads as x would underestimate the original depth variation.
# MR_cor column restricted to Observed ASVs only — residuals of Shannon,
# Simpson and Chao1 are not interpretable as index values.
#
# Reference: Schiaffino et al. (2016) Microbial Ecology 71:819-832
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 3: Library size vs richness\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
if (!requireNamespace("vegan",     quietly = TRUE)) install.packages("vegan")
library(ggplot2)
library(patchwork)
library(scales)
library(vegan)

# --- Panel builder ------------------------------------------------------------

#' Build one scatter panel: Observed ASVs or residuals vs raw read depth.
#'
#' reads and vals must be pre-matched to the same sample order.
#' Legend is hidden — shared legend added at figure level.
#'
#' @param reads     Numeric vector. Raw read depths (x-axis).
#' @param vals      Numeric vector. Observed ASVs or residuals (y-axis).
#' @param marker    Character. Marker name (drives colour).
#' @param col_title Character. Panel title.
#' @param y_label   Character. Y-axis label.
#' @param x_label   Character. X-axis label.
#' @return A ggplot object (legend hidden).
make_richness_panel <- function(reads, vals, marker,
                                col_title, y_label,
                                x_label = "Raw reads per sample") {
  
  col <- MARKER_COLORS[marker]
  df  <- data.frame(
    reads = as.numeric(reads),
    value = as.numeric(vals)
  )
  df <- df[!is.na(df$value) & !is.na(df$reads), ]
  
  if (nrow(df) < 3)
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Insufficient data",
                      color = "grey50", size = 4) +
             theme_void() +
             labs(title = col_title))
  
  # Spearman correlation
  test      <- cor.test(df$value, df$reads, method = "spearman", exact = FALSE)
  rho       <- round(test$estimate, 3)
  pval      <- round(test$p.value,  4)
  p_str     <- ifelse(pval < 0.001, "p<0.001", paste0("p=", pval))
  label     <- paste0("rho=", rho, "\n", p_str)
  ann_color <- ifelse(pval <= ALPHA_SIG, "#E84855", "grey30")
  ann_face  <- ifelse(pval <= ALPHA_SIG, "bold",    "plain")
  
  ggplot(df, aes(x = reads, y = value)) +
    geom_smooth(method = "lm", se = TRUE,
                color     = col,
                fill      = col,
                linewidth = 0.8,
                alpha     = 0.12) +
    geom_point(color = col, alpha = 0.7, size = 1.8) +
    annotate("text",
             x        = max(df$reads),
             y        = max(df$value, na.rm = TRUE),
             label    = label,
             hjust    = 1.05, vjust = 1.2,
             size     = 3,
             color    = ann_color,
             fontface = ann_face) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = col_title, x = x_label, y = y_label) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(size = 9, face = "bold", color = col),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# --- Build panels per marker --------------------------------------------------

cat("--- Building panels ----------------------------------\n\n")

all_panels <- do.call(c, lapply(MARKERS, function(m) {
  
  reads_raw <- rowSums(DS[[m]]$MR)
  
  # --- Panel 1: MR_clean — Observed ASVs vs raw reads ----------------------
  obs_clean   <- vegan::specnumber(DS[[m]]$MR_clean)
  samp_clean  <- rownames(DS[[m]]$MR_clean)
  
  p1 <- make_richness_panel(
    reads     = reads_raw[match(samp_clean, names(reads_raw))],
    vals      = as.numeric(obs_clean),
    marker    = m,
    col_title = paste0(MARKER_LABELS[m], " — Filtered (MR_clean)"),
    y_label   = "Observed ASVs"
  )
  
  # --- Panel 2: MR_rare — Observed ASVs vs raw reads -----------------------
  obs_rare   <- vegan::specnumber(DS[[m]]$MR_rare)
  samp_rare  <- rownames(DS[[m]]$MR_rare)
  
  p2 <- make_richness_panel(
    reads     = reads_raw[match(samp_rare, names(reads_raw))],
    vals      = as.numeric(obs_rare),
    marker    = m,
    col_title = paste0(MARKER_LABELS[m], " — Rarefied (MR_rare)"),
    y_label   = "Observed ASVs",
    x_label   = "Raw reads per sample (original depth)"
  )
  
  # --- Panel 3: MR_cor — log-residuals vs raw reads -------------------------
  if (!is.null(DS[[m]]$MR_cor)) {
    samp_cor  <- names(DS[[m]]$MR_cor)
    
    p3 <- make_richness_panel(
      reads     = reads_raw[match(samp_cor, names(reads_raw))],
      vals      = as.numeric(DS[[m]]$MR_cor),
      marker    = m,
      col_title = paste0(MARKER_LABELS[m], " — Corrected (MR_cor)"),
      y_label   = "Residuals (log-corrected Observed ASVs)"
    )
  } else {
    p3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "MR_cor not available",
               color = "grey50", size = 4) +
      theme_void() +
      labs(title = paste0(MARKER_LABELS[m], " — MR_cor"))
  }
  
  cat(sprintf("  %-6s | MR_clean=%d | MR_rare=%d samples\n",
              m, length(obs_clean), length(obs_rare)))
  
  list(p1, p2, p3)
}))

cat("\n")

# --- Assemble figure ----------------------------------------------------------

fig3 <- wrap_plots(all_panels,
                   nrow = length(MARKERS),
                   ncol = 3) +
  plot_annotation(
    title    = paste0(PROJECT_NAME,
                      " — Raw read depth vs Observed ASV richness"),
    subtitle = paste0(
      "x = raw reads before any filter | ",
      "Col 1: filtered | Col 2: rarefied | ",
      "Col 3: log-residuals (Observed ASVs only) | ",
      "Red rho = significant (p \u2264 ", ALPHA_SIG, ")"),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8,  color = "grey40"))
  )

print(fig3)
save_figure(fig3, "02_03_library_size_vs_richness",
            width = 14, height = 4 * length(MARKERS) + 1)

cat("  Figure 3 saved to:", PATH_FIG_DIR, "\n\n")


# =============================================================================
# SECTION 4 — RAREFACTION CURVES
# =============================================================================
# Rarefaction curves showing ASV richness as a function of sequencing depth
# for all samples. One figure per marker.
# Based on MR_clean (filtered matrix) — the input to rarefaction.
# Rarefaction threshold overlaid as vertical dashed line.
#
# Purpose: assess whether sequencing depth is sufficient to capture community
# richness and validate the chosen rarefaction threshold.
#
# Reference: Oksanen et al. (2024) vegan. R package.
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 4: Rarefaction curves\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("vegan",     quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
library(vegan)
library(ggplot2)
library(scales)

# --- Build rarefaction curve data frame per marker ----------------------------

#' Compute rarefaction curves for one marker and return a long data.frame.
#'
#' Uses vegan::rarecurve() internally and reformats output for ggplot.
#'
#' @param mr     data.frame. Samples x ASVs count matrix (MR_clean).
#' @param marker Character. Marker name.
#' @param step   Integer. Step size for rarefaction interpolation.
#' @return Long-format data.frame with columns: sample, reads, richness.
build_rarecurve_df <- function(mr, marker, step = 100) {
  
  cat(sprintf("  Computing rarefaction curves for %s (%d samples)...\n",
              marker, nrow(mr)))
  
  # tmax set to 3x rarefaction threshold — covers the informative range
  # without computing up to max read depth (slow for deeply sequenced samples)
  tmax_val <- RARE_THRESHOLDS[[marker]] * 3
  rc <- vegan::rarecurve(mr, step = step,
                         tmax  = tmax_val,
                         label = FALSE)
  
  do.call(rbind, lapply(seq_along(rc), function(i) {
    samp  <- rownames(mr)[i]
    reads <- attr(rc[[i]], "Subsample")
    data.frame(
      sample     = samp,
      reads      = as.numeric(reads),
      richness   = as.numeric(rc[[i]]),
      stringsAsFactors = FALSE
    )
  }))
}

# --- Build figures per marker -------------------------------------------------

for (m in MARKERS) {
  
  mr     <- DS[[m]]$MR_clean
  thresh <- RARE_THRESHOLDS[[m]]
  col    <- MARKER_COLORS[m]
  step   <- RARE_STEP[[m]]
  
  # Compute curves
  rc_df <- build_rarecurve_df(mr, marker = m, step = step)

  # Samples below threshold
  max_reads     <- tapply(rc_df$reads, rc_df$sample, max)
  n_below       <- sum(max_reads < thresh)
  n_total       <- nrow(mr)
  
  subtitle <- sprintf(
    "%d samples | Threshold = %s reads | %d sample(s) below threshold",
    n_total,
    format(thresh, big.mark = ","),
    n_below)
  
  fig <- ggplot(rc_df,
                aes(x     = reads,
                    y     = richness,
                    group = sample)) +
    geom_line(color = col, alpha = 0.6, linewidth = 0.5) +
    geom_vline(xintercept = thresh,
               linetype   = "dashed",
               color      = "grey20",
               linewidth  = 0.8) +
    annotate("text",
             x     = thresh,
             y     = Inf,
             label = paste0("threshold\n",
                            format(thresh, big.mark = ",")),
             hjust = -0.1, vjust = 1.3,
             size  = 3, color = "grey20") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = paste0(MARKER_LABELS[m],
                        " — Rarefaction curves"),
      subtitle = subtitle,
      x        = "Reads per sample",
      y        = "Observed ASVs"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(size = 11, face = "bold", color = col),
      plot.subtitle    = element_text(size = 9,  color = "grey40"),
      panel.grid.minor = element_blank()
    )
  
  print(fig)
  save_figure(fig,
              paste0("02_04_rarefaction_curves_", m),
              width = 10, height = 7)
  cat(sprintf("  Figure 4 (%s) saved.\n\n", m))
}

# =============================================================================
# SECTION 5 — ASV SEQUENCE LENGTH DISTRIBUTION
# =============================================================================
# Histograms of ASV sequence lengths from the filtered community (MR_clean +
# ASS_clean). Two weighting approaches shown per marker:
#   - ASV-weighted  : each ASV counts once regardless of read abundance
#   - Read-weighted : each ASV weighted by its total read count
#
# Expected amplicon length range overlaid as shaded region (AMPLICON_LENGTHS).
# ASVs outside the expected range may represent non-target amplicons or
# chimeras not removed by DADA2.
#
# One figure per marker (two panels: ASV-weighted and read-weighted).
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 5: ASV sequence length distribution\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
library(ggplot2)
library(patchwork)
library(scales)

# --- Check for sequence length column in ASS ----------------------------------
# ASV sequence lengths must be stored in ASS_clean.
# Column name checked — common names: "Length", "length", "seq_length".

# Check for Sequence column — compute lengths on the fly
SEQ_COL <- NULL
for (col_candidate in c("Sequence", "sequence", "seq")) {
  if (all(sapply(MARKERS, function(m)
    col_candidate %in% colnames(DS[[m]]$ASS_clean)))) {
    SEQ_COL <- col_candidate
    break
  }
}

if (is.null(SEQ_COL)) {
  cat("  [SKIP] No sequence column found in ASS_clean.\n")
  cat("  Checked: Sequence, sequence, seq\n")
  cat("  Section 5 skipped.\n\n")
} else {
  cat(sprintf("  Sequence column found: '%s' — computing lengths\n\n", SEQ_COL))
  
  # --- Build length data frame ------------------------------------------------
  
  build_length_df <- function(marker) {
    
    ass   <- DS[[marker]]$ASS_clean
    mr    <- DS[[marker]]$MR_clean
    
    # Total reads per ASV across all samples
    asv_reads <- colSums(mr)
    
    # Match ASVs between ASS and MR — should be identical after syncing
    common_asvs <- intersect(rownames(ass), names(asv_reads))
    
    data.frame(
      asv    = common_asvs,
      marker = marker,
      length = nchar(as.character(ass[common_asvs, SEQ_COL])),
      reads  = as.numeric(asv_reads[common_asvs]),
      stringsAsFactors = FALSE
    )
  }
  
  length_list <- setNames(lapply(MARKERS, build_length_df), MARKERS)
  
  # Console summary
  cat("--- Length distribution summary ----------------------\n\n")
  for (m in MARKERS) {
    df     <- length_list[[m]]
    lrange <- AMPLICON_LENGTHS[[m]]
    n_out  <- if (!is.null(lrange))
      sum(df$length < lrange[1] | df$length > lrange[2], na.rm = TRUE)
    else NA
    
    cat(sprintf("  %-6s | %d ASVs | length: %d-%d bp | median=%d bp",
                m, nrow(df),
                min(df$length, na.rm = TRUE),
                max(df$length, na.rm = TRUE),
                median(df$length, na.rm = TRUE)))
    if (!is.na(n_out))
      cat(sprintf(" | %d ASVs outside expected range [%d-%d bp]",
                  n_out, lrange[1], lrange[2]))
    cat("\n")
  }
  cat("\n")
  
  # --- Build figure per marker ------------------------------------------------
  
  for (m in MARKERS) {
    
    df     <- length_list[[m]]
    col    <- MARKER_COLORS[m]
    lrange <- AMPLICON_LENGTHS[[m]]
    
    # --- Panel A: ASV-weighted ------------------------------------------------
    pA <- ggplot(df, aes(x = length)) +
      geom_histogram(fill = col, color = "white",
                     binwidth = 5, alpha = 0.8) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma) +
      labs(title    = "ASV-weighted",
           subtitle = "Each ASV counted once",
           x        = "Sequence length (bp)",
           y        = "Number of ASVs") +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(size = 10, face = "bold"),
            plot.subtitle    = element_text(size = 8,  color = "grey40"),
            panel.grid.minor = element_blank())
    
    # --- Panel B: Read-weighted -----------------------------------------------
    pB <- ggplot(df, aes(x = length, weight = reads)) +
      geom_histogram(fill = col, color = "white",
                     binwidth = 5, alpha = 0.8) +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma) +
      labs(title    = "Read-weighted",
           subtitle = "Each ASV weighted by total read count",
           x        = "Sequence length (bp)",
           y        = "Number of reads") +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(size = 10, face = "bold"),
            plot.subtitle    = element_text(size = 8,  color = "grey40"),
            panel.grid.minor = element_blank())
    
    # Add expected length range shading if defined
    if (!is.null(lrange)) {
      shade <- annotate("rect",
                        xmin  = lrange[1], xmax = lrange[2],
                        ymin  = -Inf,      ymax = Inf,
                        fill  = "grey80",  alpha = 0.3)
      range_label <- annotate("text",
                              x     = mean(lrange),
                              y     = Inf,
                              label = paste0("Expected range\n",
                                             lrange[1], "–", lrange[2], " bp"),
                              vjust = 1.3, size = 3, color = "grey30")
      pA <- pA + shade + range_label
      pB <- pB + shade + range_label
    }
    
    fig5 <- (pA | pB) +
      plot_annotation(
        title    = paste0(MARKER_LABELS[m],
                          " — ASV sequence length distribution (MR_clean)"),
        subtitle = sprintf("%d ASVs | %s total reads",
                           nrow(df),
                           format(sum(df$reads), big.mark = ",")),
        theme = theme(
          plot.title    = element_text(size = 12, face = "bold", color = col),
          plot.subtitle = element_text(size = 9,  color = "grey40"))
      )
    
    print(fig5)
    save_figure(fig5,
                paste0("02_05_asv_length_", m),
                width = 12, height = 5)
    cat(sprintf("  Figure 5 (%s) saved.\n", m))
  }
  cat("\n")
  
} # end SEQ_COL check

# =============================================================================
# SECTION 6 — IDENTITY DISTRIBUTION
# =============================================================================
# Histograms of VSEARCH % identity values from the filtered community
# (ASS_clean). Shows the distribution of taxonomic assignment quality.
#
# Two panels per marker:
#   Left  : raw counts (number of ASVs per identity bin)
#   Right : read-weighted (number of reads per identity bin)
#
# x-axis starts at filter threshold — ASVs below threshold already excluded
# by 01_filtering.R. The axis start implicitly documents the threshold value.
#
# Approach after Mahé et al. (2017) — identity distribution used to assess
# the reliability of taxonomic assignments and validate filter thresholds.
#
# Reference: Mahé et al. (2017) Swarm v2. PeerJ 5:e3166
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 6: Identity distribution\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
library(ggplot2)
library(patchwork)
library(scales)

# --- Parameters ---------------------------------------------------------------
# x-axis starts at filter threshold — ASVs below threshold already excluded
# by 01_filtering.R. The axis start implicitly documents the threshold value.

#' Retrieve identity filter threshold for one marker from DS$PARAMS.
#'
#' @param marker Character. Marker name.
#' @return Numeric. Identity threshold percentage (default 60 if not found).
get_identity_threshold <- function(marker) {
  t <- DS[[marker]]$PARAMS$filter_min_identity
  if (is.null(t)) {
    warning(marker, ": PARAMS$filter_min_identity not found — using default 60")
    return(60)
  }
  return(t)
}

# --- Build identity data frame ------------------------------------------------

#' Build identity distribution data frame for one marker.
#'
#' @param marker Character. Marker name.
#' @return A data.frame with columns: asv, identity, reads.
build_identity_df <- function(marker) {
  
  ass <- DS[[marker]]$ASS_clean
  mr  <- DS[[marker]]$MR_clean
  
  # Total reads per ASV
  asv_reads   <- colSums(mr)
  common_asvs <- intersect(rownames(ass), names(asv_reads))
  
  if (!"Pct_identity" %in% colnames(ass))
    stop(marker, ": 'Pct_identity' column not found in ASS_clean")
  
  data.frame(
    asv      = common_asvs,
    marker   = marker,
    identity = as.numeric(ass[common_asvs, "Pct_identity"]),
    reads    = as.numeric(asv_reads[common_asvs]),
    stringsAsFactors = FALSE
  )
}

identity_list <- setNames(lapply(MARKERS, build_identity_df), MARKERS)

# --- Console summary ----------------------------------------------------------

cat("--- Identity distribution summary -------------------\n\n")
for (m in MARKERS) {
  df     <- identity_list[[m]]
  thresh <- get_identity_threshold(m)
  n_out  <- sum(df$identity < thresh, na.rm = TRUE)
  cat(sprintf(
    "  %-6s | %d ASVs | identity: %.1f%%-%.1f%% | median=%.1f%% | <%.0f%%: %d ASVs\n",
    m, nrow(df),
    min(df$identity,    na.rm = TRUE),
    max(df$identity,    na.rm = TRUE),
    median(df$identity, na.rm = TRUE),
    thresh, n_out))
}
cat("\n")

# --- Build figure per marker --------------------------------------------------

for (m in MARKERS) {
  
  df     <- identity_list[[m]]
  col    <- MARKER_COLORS[m]
  thresh <- get_identity_threshold(m)

  # --- Panel A: ASV-weighted --------------------------------------------------
  pA <- ggplot(df, aes(x = identity)) +
    geom_histogram(fill     = col,
                   color    = "white",
                   binwidth = IDENTITY_BINWIDTH,
                   alpha    = 0.8) +
    scale_x_continuous(limits = c(thresh, 100),
                       oob    = scales::squish,
                       labels = function(x) paste0(x, "%")) +
    scale_y_continuous(labels = scales::comma) +
    labs(title    = "ASV-weighted",
         subtitle = "Each ASV counted once",
         x        = "% identity to reference",
         y        = "Number of ASVs") +
    theme_bw(base_size = 10) +
    theme(plot.title       = element_text(size = 10, face = "bold"),
          plot.subtitle    = element_text(size = 8,  color = "grey40"),
          panel.grid.minor = element_blank())
  
  # --- Panel B: Read-weighted -------------------------------------------------
  pB <- ggplot(df, aes(x = identity, weight = reads)) +
    geom_histogram(fill     = col,
                   color    = "white",
                   binwidth = IDENTITY_BINWIDTH,
                   alpha    = 0.8) +
    scale_x_continuous(limits = c(thresh, 100),
                       oob    = scales::squish,
                       labels = function(x) paste0(x, "%")) +
    scale_y_continuous(labels = scales::comma) +
    labs(title    = "Read-weighted",
         subtitle = "Each ASV weighted by total read count",
         x        = "% identity to reference",
         y        = "Number of reads") +
    theme_bw(base_size = 10) +
    theme(plot.title       = element_text(size = 10, face = "bold"),
          plot.subtitle    = element_text(size = 8,  color = "grey40"),
          panel.grid.minor = element_blank())
  
  fig6 <- (pA | pB) +
    plot_annotation(
      title    = paste0(MARKER_LABELS[m],
                        " — % identity distribution (MR_clean)"),
      subtitle = sprintf(
        "%d ASVs | %s total reads | x-axis starts at filter threshold (%.0f%%)",
        nrow(df),
        format(sum(df$reads), big.mark = ","),
        thresh),
      theme = theme(
        plot.title    = element_text(size = 12, face = "bold", color = col),
        plot.subtitle = element_text(size = 9,  color = "grey40"))
    )
  
  print(fig6)
  save_figure(fig6,
              paste0("02_06_identity_distribution_", m),
              width = 12, height = 5)
  cat(sprintf("  Figure 6 (%s) saved.\n", m))
}
cat("\n")

# =============================================================================
# SECTION 7 — FILTERING JOURNEY VISUALISATION
# =============================================================================
# Documents the step-by-step impact of each filtering stage on sample count,
# ASV richness and read depth. One table and one figure per marker.
#
# Steps tracked:
#   1. Raw          : $MR (all samples, all ASVs, before any filter)
#   2. Filtered     : $MR_clean (after taxonomic + rare ASV + quality + exclusions)
#   3. Rarefied     : $MR_rare (after rarefaction threshold)
#
# For 18S, the taxonomic filter (non-protist removal) is the dominant step
# and is documented separately via DS$PARAMS if available.
#
# Output:
#   CSV  : one row per marker x step with n_samples, n_asvs, n_reads, % retained
#   Fig  : horizontal bar chart of % retained at each step per marker
# =============================================================================

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 7: Filtering journey\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
library(ggplot2)
library(patchwork)
library(scales)

# --- Build filtering journey table --------------------------------------------

build_journey <- function(marker) {
  
  # Raw values
  mr_raw    <- DS[[marker]]$MR
  mr_clean  <- DS[[marker]]$MR_clean
  mr_rare   <- DS[[marker]]$MR_rare
  
  n_samp_raw   <- nrow(mr_raw)
  n_asvs_raw   <- ncol(mr_raw)
  n_reads_raw  <- sum(mr_raw)
  
  n_samp_clean  <- nrow(mr_clean)
  n_asvs_clean  <- ncol(mr_clean)
  n_reads_clean <- sum(mr_clean)
  
  n_samp_rare   <- nrow(mr_rare)
  n_asvs_rare   <- ncol(mr_rare)
  n_reads_rare  <- sum(mr_rare)
  
  data.frame(
    marker   = marker,
    step     = c("Raw", "Filtered", "Rarefied"),
    n_samples = c(n_samp_raw,   n_samp_clean,  n_samp_rare),
    n_asvs    = c(n_asvs_raw,   n_asvs_clean,  n_asvs_rare),
    n_reads   = c(n_reads_raw,  n_reads_clean, n_reads_rare),
    pct_samp  = round(c(100,
                        n_samp_clean  / n_samp_raw  * 100,
                        n_samp_rare   / n_samp_raw  * 100), 1),
    pct_asvs  = round(c(100,
                        n_asvs_clean  / n_asvs_raw  * 100,
                        n_asvs_rare   / n_asvs_raw  * 100), 1),
    pct_reads = round(c(100,
                        n_reads_clean / n_reads_raw  * 100,
                        n_reads_rare  / n_reads_raw  * 100), 1),
    stringsAsFactors = FALSE
  )
}

journey_all <- do.call(rbind, lapply(MARKERS, build_journey))
journey_all$marker <- factor(journey_all$marker, levels = MARKERS)
journey_all$step   <- factor(journey_all$step,
                             levels = c("Raw", "Filtered", "Rarefied"))

# --- Console summary ----------------------------------------------------------

cat("--- Filtering journey summary ------------------------\n\n")
cat(sprintf("  %-6s | %-10s | %8s | %8s | %12s | %8s | %8s | %8s\n",
            "Marker", "Step", "Samples", "ASVs", "Reads",
            "%Samp", "%ASVs", "%Reads"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (m in MARKERS) {
  sub <- journey_all[journey_all$marker == m, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("  %-6s | %-10s | %8d | %8d | %12s | %7.1f%% | %7.1f%% | %7.1f%%\n",
                r$marker, r$step,
                r$n_samples, r$n_asvs,
                format(r$n_reads, big.mark = ","),
                r$pct_samp, r$pct_asvs, r$pct_reads))
  }
  cat(paste(rep("-", 80), collapse = ""), "\n")
}
cat("\n")

# Save CSV
PATH_JOURNEY_CSV <- file.path(PATH_TAB_DIR,
                              paste0(PROJECT_NAME,
                                     "_02_filtering_journey.csv"))
write.table(journey_all, file = PATH_JOURNEY_CSV,
            sep = ";", row.names = FALSE, quote = TRUE)
cat("  Journey table saved to:", PATH_JOURNEY_CSV, "\n\n")

# --- Figure 7 — % retained at each step per marker ---------------------------

cat("--- Figure 7: Filtering journey ---------------------\n")

# Reshape to long format for plotting
metric_labels <- c(pct_samp = "Samples", pct_asvs = "ASVs", pct_reads = "Reads")

plot_df <- do.call(rbind, lapply(names(metric_labels), function(metric) {
  data.frame(
    marker = journey_all$marker,
    step   = journey_all$step,
    metric = metric_labels[[metric]],
    pct    = journey_all[[metric]],
    stringsAsFactors = FALSE
  )
}))

plot_df$metric <- factor(plot_df$metric,
                         levels = c("Samples", "ASVs", "Reads"))
plot_df$step   <- factor(plot_df$step,
                         levels = c("Raw", "Filtered", "Rarefied"))

fig7_list <- lapply(MARKERS, function(m) {
  
  sub <- plot_df[plot_df$marker == m, ]
  col <- MARKER_COLORS[m]
  
  ggplot(sub, aes(x = step, y = pct, fill = metric)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.6, alpha = 0.85) +
    geom_text(aes(label = paste0(pct, "%")),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 2.8, color = "grey20") +
    scale_fill_manual(
      values = c("Samples" = "#4E79A7",
                 "ASVs"    = "#59A14F",
                 "Reads"   = "#F28E2B"),
      name = NULL) +
    scale_y_continuous(limits = c(0, 115),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = MARKER_LABELS[m],
         x     = NULL,
         y     = "% retained vs raw") +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(size = 10, face = "bold", color = col),
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
})

fig7 <- wrap_plots(fig7_list, nrow = 1) +
  plot_annotation(
    title    = paste0(PROJECT_NAME, " — Filtering journey"),
    subtitle = paste0(
      "% retained relative to raw at each pipeline stage | ",
      "Raw = 100% by definition"),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40"))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(fig7)
save_figure(fig7, "02_07_filtering_journey",
            width = 14, height = 6)
cat("  Figure 7 saved.\n\n")

# =============================================================================
# SECTION 8 — SESSION INFO
# =============================================================================
# Saves R session information for reproducibility reporting.

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 8: Session info\n")
cat("  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")

sink(PATH_SESSION_TXT)
cat(PROJECT_NAME, "— 02_dataset_quality.R session information\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

cat("  Saved:", PATH_SESSION_TXT, "\n\n")


# =============================================================================
# SECTION 9 — CLEAN ENVIRONMENT
# =============================================================================
# Retains only DS and project identifiers for downstream scripts.
# All intermediate objects, figures, and working variables are removed.

cat("=============================================================================\n")
cat(" ", PROJECT_NAME, "— Section 9: Clean environment\n")
cat("=============================================================================\n\n")

rm(list = setdiff(ls(), c(
  "DS",
  "ROOT_DIR",
  "PROJECT_NAME",
  "MARKERS",
  "MARKER_LABELS",
  "MARKER_COLORS"
)))

cat("  Environment cleaned.\n")
cat("  Objects retained: DS, ROOT_DIR, PROJECT_NAME,\n")
cat("                    MARKERS, MARKER_LABELS, MARKER_COLORS\n")
cat("\n  Statistics stage complete — diagnostic figures and tables written.\n")
cat("=============================================================================\n")

