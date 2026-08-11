# =============================================================================
# File      : config.R
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Marker    : 18S (Protists)
# Author    : David Singer
# Date      : 2026-03-03
# Description: DADA2 and taxonomy parameters for 18S marker.
#              Sourced by 04_DADA2_ASV.R and 05_taxonomy.R.
#              PROJECT_ROOT is derived from this file's path in the R scripts.
# =============================================================================

# --- Marker ------------------------------------------------------------------
MARKER            <- "18S"
TAXONOMIC_GROUP   <- "protists"

# --- Directories -------------------------------------------------------------
# config_path is passed by 04_DADA2_ASV.R / 05_taxonomy.R before source()
.project_root     <- normalizePath(file.path(dirname(config_path), ".."))
PIPELINE_DIR      <- file.path(.project_root, "results", MARKER)
STEP3_DIR         <- file.path(PIPELINE_DIR, "03_trimmed")
STEP4_DIR         <- file.path(PIPELINE_DIR, "04_ASV")
STEP5_DIR         <- file.path(PIPELINE_DIR, "05_taxonomy")
LOG_DIR           <- file.path(PIPELINE_DIR, "logs")
FILTERED_DIR      <- file.path(STEP4_DIR,    "tmp_filtered")
READ_TRACKING     <- file.path(LOG_DIR,       "read_tracking.tsv")

# --- Databases ---------------------------------------------------------------
DB_PR2_VSEARCH    <- file.path(.project_root, "databases", "PR2", "pr2_version_5.1.0_SSU_mothur.fasta")
DB_PR2_DADA2      <- file.path(.project_root, "databases", "PR2", "pr2_version_5.1.0_SSU_dada2.fasta")

# --- File naming -------------------------------------------------------------
FASTQ_STRIP_PATTERN_R <- "_R1\\.fastq\\.gz$"   # regex to strip R1 suffix → sample name

# --- DADA2 filterAndTrim parameters ------------------------------------------
TRUNCLEN_FWD      <- 275
TRUNCLEN_REV      <- 275
MAX_EE_FWD        <- 1
MAX_EE_REV        <- 1
TRUNC_Q           <- 2
MAX_N             <- 0
COMPRESS          <- TRUE

# --- DADA2 denoising parameters ----------------------------------------------
POOLING <- FALSE     # "pseudo", "pooled", or FALSE
CHIMERA_METHOD    <- "consensus"  # "consensus", "pooled", or "per-sample"
THREADS           <- 32
THREADS_LEARN     <- 1            # learnErrors only — keep at 1 to avoid Linux segfault
MIN_OVERLAP       <- 8            # 18S V4 ~380bp, variable amplicon — relax overlap

# --- VSEARCH parameters ------------------------------------------------------
VSEARCH_PATH      <- "vsearch"    # full path if not in PATH
VSEARCH_ID        <- 0.1
VSEARCH_MAXACCEPTS <- 1
VSEARCH_MAXREJECTS <- 32
VSEARCH_THREADS   <- 32

# --- Taxonomy levels ---------------------------------------------------------
# PR2 v5.1.0 has 9 levels — assignTaxonomy() returns generic names,
# we rename manually using TAXO_LEVELS after assignment.
# Unified 9-level standard: Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species
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

# PR2 v5.1.0 mothur header format: >AB353770.1.1740_U
# GenBank extraction: strip trailing length suffix e.g. ".1740_U"
# Handled in 05_taxonomy.R

# --- Output filenames --------------------------------------------------------
OUT_FASTA         <- file.path(STEP4_DIR, paste0(MARKER, "_Fasta.fasta"))
OUT_MR            <- file.path(STEP4_DIR, paste0(MARKER, "_MR.csv"))
OUT_ASS           <- file.path(STEP5_DIR, paste0(MARKER, "_ASS.csv"))
OUT_SEQTAB_RDS    <- file.path(STEP4_DIR, "seqtab_nochim.rds")
OUT_VSEARCH_BLAST6 <- file.path(STEP5_DIR, paste0(MARKER, "_vsearch.blast6"))
