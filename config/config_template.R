# =============================================================================
# File      : config_template.R
# Author    : David Singer
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Description: Template for marker-specific R configuration.
#              Copy this file to your Pipeline_<MARKER>/ directory, rename it
#              config.R, and fill in all PLACEHOLDER values.
#              Sourced by 04_DADA2_ASV.R and 05_taxonomy.R.
# Usage     : cp config/config_template.R Pipeline_MYMARKER/config.R
#             # then edit the copy — do NOT edit this template directly
# =============================================================================

# =============================================================================
# SECTION 1 — Marker identity
# =============================================================================

# Short marker label — must match MARKER in config.sh for the same pipeline
MARKER           <- "MARKER_PLACEHOLDER"   # e.g. "18S", "16S", "ITS", "COI"

# Free-text description of the target group — used in log/report headers only
TAXONOMIC_GROUP  <- "your_taxonomic_group" # e.g. "protists", "bacteria", "fungi"

# =============================================================================
# SECTION 2 — Directories (do not modify — derived from config_path)
# =============================================================================
# config_path is set by 04_DADA2_ASV.R and 05_taxonomy.R before calling
# source(config_path), making this file's own location available at parse time.

.project_root  <- normalizePath(file.path(dirname(config_path), ".."))
PIPELINE_DIR   <- file.path(.project_root, "results", MARKER)
STEP3_DIR      <- file.path(PIPELINE_DIR, "03_trimmed")
STEP4_DIR      <- file.path(PIPELINE_DIR, "04_ASV")
STEP5_DIR      <- file.path(PIPELINE_DIR, "05_taxonomy")
LOG_DIR        <- file.path(PIPELINE_DIR, "logs")
FILTERED_DIR   <- file.path(STEP4_DIR,   "tmp_filtered")
READ_TRACKING  <- file.path(LOG_DIR,      "read_tracking.tsv")

# =============================================================================
# SECTION 3 — Reference databases
# =============================================================================
# Provide paths relative to project root using the helper below.
# Never hardcode absolute paths — this keeps the config portable.
#
# Example database entries by marker:
#
#   18S (PR2 v5.1.0):
#     DB_VSEARCH <- file.path(.project_root, "databases", "PR2",
#                             "pr2_version_5.1.0_SSU_mothur.fasta")
#     DB_DADA2   <- file.path(.project_root, "databases", "PR2",
#                             "pr2_version_5.1.0_SSU_dada2.fasta")
#
#   ITS (UNITE v9 / 19.02.2025):
#     DB_VSEARCH <- file.path(.project_root, "databases", "UNITE",
#                             "sh_general_release_dynamic_19.02.2025.fasta")
#     DB_DADA2   <- file.path(.project_root, "databases", "UNITE",
#                             "sh_general_release_dynamic_19.02.2025_dada2_fixed.fasta")
#
#   16S (SILVA 138.2):
#     DB_VSEARCH <- file.path(.project_root, "databases", "SILVA",
#                             "SILVA_138.2_SSURef_NR99_tax_silva.fasta")
#     DB_DADA2   <- file.path(.project_root, "databases", "SILVA",
#                             "silva_nr99_v138.2_train_set.fa.gz")

DB_VSEARCH     <- file.path(.project_root, "databases", "YOUR_DB", "vsearch_database.fasta")
DB_DADA2       <- file.path(.project_root, "databases", "YOUR_DB", "dada2_database.fasta")

# =============================================================================
# SECTION 4 — File naming
# =============================================================================

# Regex used to strip the R1 suffix from FASTQ filenames to derive sample names.
# Matches the end of the filename, e.g. "Sample1_R1.fastq.gz" -> "Sample1".
# Adjust the pattern if your filenames use a different convention.
FASTQ_STRIP_PATTERN_R <- "_R1\\.fastq\\.gz$"

# =============================================================================
# SECTION 5 — DADA2: filterAndTrim parameters
# =============================================================================
# Review trimmed read quality profiles (MultiQC step 2 output) before setting
# truncation lengths. Aim to retain the highest-quality portion of each read
# while still allowing forward and reverse reads to overlap after merging.

# Truncation position for forward reads (bp). Set to 0 to disable truncation.
# Disable (0) for variable-length amplicons such as ITS.
TRUNCLEN_FWD   <- 230

# Truncation position for reverse reads (bp). Set to 0 to disable truncation.
# Reverse reads often degrade sooner — a shorter value is typical.
TRUNCLEN_REV   <- 180

# Maximum expected errors (EE) allowed per read after truncation.
# Reads exceeding this threshold are discarded. Lower = stricter.
# Typical values: 1–2 for high-quality runs; relax to 3–5 for older data.
MAX_EE_FWD     <- 2
MAX_EE_REV     <- 2

# Truncate at first base below this quality score (Phred).
# 2 = DADA2 default — removes only the worst-quality tails.
TRUNC_Q        <- 2

# Maximum number of ambiguous (N) bases allowed per read. 0 = reject any N.
MAX_N          <- 0

# Compress filtered output files (TRUE = gzip). Recommended to save disk space.
COMPRESS       <- TRUE

# =============================================================================
# SECTION 6 — DADA2: denoising parameters
# =============================================================================

# Sample pooling strategy for dada():
#   FALSE    = independent per-sample denoising (fastest, default)
#   "pseudo" = pseudo-pooling — improves rare variant detection, moderate cost
#   "pooled" = full pooling — best sensitivity, highest memory/time cost
POOLING        <- FALSE

# Chimera detection method passed to removeBimeraDenovo():
#   "consensus"  = majority-vote per-sample (default, fast)
#   "pooled"     = across all samples (slower, slightly more sensitive)
#   "per-sample" = per sample (equivalent to "consensus" for single runs)
CHIMERA_METHOD <- "consensus"

# Number of threads for dada() and mergePairs(). Set to 1 on Windows.
THREADS        <- 4

# Number of threads for learnErrors() ONLY.
# Keep at 1 on Linux to avoid a known multithreaded segfault in DADA2 < 1.28.
THREADS_LEARN  <- 1

# Minimum overlap (bp) required for forward/reverse read merging.
# Must be ≥ TRUNCLEN_FWD + TRUNCLEN_REV − amplicon_length.
# Typical values:
#   18S V4  (~380 bp amplicon, 275+275 reads): 8–12
#   ITS2    (150–350 bp variable amplicon)   : 8
#   16S V4-V5 (~411 bp, 275+275 reads)       : 12–20
MIN_OVERLAP    <- 12

# =============================================================================
# SECTION 7 — VSEARCH parameters
# =============================================================================

# Path to vsearch executable. Use "vsearch" if it is on PATH, otherwise
# provide the full absolute path (e.g. "/home/user/miniconda3/bin/vsearch").
VSEARCH_PATH      <- "vsearch"

# Minimum sequence identity for a VSEARCH hit to be accepted (0.0–1.0).
# Note: VSEARCH uses global alignment — 0.1 accepts very distant matches
# to ensure all ASVs get assigned; post-filter by confidence/score instead.
VSEARCH_ID        <- 0.1

# Maximum number of accepted hits per query. 1 returns the best hit only.
VSEARCH_MAXACCEPTS <- 1

# Maximum number of rejected candidates before giving up on a query.
# Higher values increase sensitivity at the cost of speed. 32 is standard.
VSEARCH_MAXREJECTS <- 32

# Number of threads for VSEARCH.
VSEARCH_THREADS   <- 4

# =============================================================================
# SECTION 8 — Taxonomy levels
# =============================================================================
# MarkerFlow uses a unified 9-level taxonomy regardless of the database.
# Levels not present in a given database are filled with NA automatically.
#
# Database-to-level mapping:
#   PR2 v5.1.0  : all 9 levels present natively
#   UNITE       : 7 levels; Domain and Supergroup hardcoded via TAXO_FIXED
#   SILVA 138.2 : 6 levels (Kingdom, Phylum, Class, Order, Family, Genus)
#                 mapped via TAXO_SILVA_MAP; Division/Subdivision/Species = NA

TAXO_LEVELS <- c(
    "Domain",
    "Supergroup",
    "Division",
    "Subdivision",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
)

# --- Optional: hardcode top-level taxa (use for single-kingdom databases) ----
# Set to NULL if your database covers multiple kingdoms (e.g. SILVA, PR2).
# Set to a named list to override fixed levels (e.g. UNITE covers only Fungi).
#
# Example for UNITE (ITS):
#   TAXO_FIXED <- list(Domain = "Eukaryota", Supergroup = "Obazoa", Division = NA)
#
# Example for PR2 / SILVA (multi-kingdom):
#   TAXO_FIXED <- NULL

TAXO_FIXED <- NULL

# --- Optional: strip database-specific prefixes from taxonomy strings --------
# Set TRUE for UNITE (e.g. "p__Ascomycota" → "Ascomycota").
# Set FALSE for PR2 and SILVA (no prefixes).
TAXO_STRIP_PREFIX <- FALSE

# --- Optional: SILVA-specific level mapping ----------------------------------
# Only needed when using a SILVA-formatted database that returns 6 levels
# named differently from the MarkerFlow 9-level standard.
# Set to NULL for non-SILVA databases.
#
# Example for SILVA 138.2:
#   TAXO_SILVA_LEVELS <- c("Domain", "Supergroup", "Class", "Order", "Family", "Genus")
#   TAXO_SILVA_MAP    <- c(Domain="Domain", Supergroup="Supergroup",
#                          Division=NA, Subdivision=NA,
#                          Class="Class", Order="Order",
#                          Family="Family", Genus="Genus", Species=NA)

TAXO_SILVA_LEVELS <- NULL
TAXO_SILVA_MAP    <- NULL

# =============================================================================
# SECTION 9 — Output filenames (do not modify)
# =============================================================================

OUT_FASTA          <- file.path(STEP4_DIR, paste0(MARKER, "_Fasta.fasta"))
OUT_MR             <- file.path(STEP4_DIR, paste0(MARKER, "_MR.csv"))
OUT_ASS            <- file.path(STEP5_DIR, paste0(MARKER, "_ASS.csv"))
OUT_SEQTAB_RDS     <- file.path(STEP4_DIR, "seqtab_nochim.rds")
OUT_VSEARCH_BLAST6 <- file.path(STEP5_DIR, paste0(MARKER, "_vsearch.blast6"))
