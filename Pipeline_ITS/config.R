# =============================================================================
# File      : config.R
# Project   : MetaBarFlow – generalised metabarcoding pipeline
# Marker    : ITS (Fungi)
# Author    : David Singer
# Date      : 2026-03-03
# Description: DADA2 and taxonomy parameters for ITS marker.
#              Sourced by 04_DADA2_ASV.R and 05_taxonomy.R.
#              PROJECT_ROOT is derived from this file's path in the R scripts.
# =============================================================================

# --- Marker ------------------------------------------------------------------
MARKER            <- "ITS"
TAXONOMIC_GROUP   <- "fungi"

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
DB_UNITE_VSEARCH  <- file.path(.project_root, "databases", "UNITE", "sh_general_release_dynamic_19.02.2025.fasta")
DB_UNITE_DADA2    <- file.path(.project_root, "databases", "UNITE", "sh_general_release_dynamic_19.02.2025_dada2_fixed.fasta")

# --- File naming -------------------------------------------------------------
FASTQ_STRIP_PATTERN_R <- "_R1\\.fastq\\.gz$"   # regex to strip R1 suffix → sample name

# --- DADA2 filterAndTrim parameters ------------------------------------------
# ITS is variable length — do NOT truncate (truncLen = 0 disables truncation)
TRUNCLEN_FWD      <- 0
TRUNCLEN_REV      <- 0
MAX_EE_FWD        <- 1
MAX_EE_REV        <- 1
TRUNC_Q           <- 2
MAX_N             <- 0
COMPRESS          <- TRUE

# --- DADA2 denoising parameters ----------------------------------------------
POOLING           <- FALSE
CHIMERA_METHOD    <- "consensus"
THREADS           <- 18
THREADS_LEARN     <- 1            # learnErrors only — keep at 1 to avoid Linux segfault
MIN_OVERLAP       <- 8            # ITS variable length 150-350bp — relax overlap

# --- VSEARCH parameters ------------------------------------------------------
VSEARCH_PATH      <- "vsearch"    # full path if not in PATH
VSEARCH_ID        <- 0.1
VSEARCH_MAXACCEPTS <- 1
VSEARCH_MAXREJECTS <- 32
VSEARCH_THREADS   <- 32

# --- Taxonomy levels ---------------------------------------------------------
# UNITE has 7 native levels — mapped to 9-level standard:
# Domain and Supergroup are hardcoded (see TAXO_FIXED below)
# Division = Phylum (p__), Subdivision = NA, Class = c__, Order = o__,
# Family = f__, Genus = g__, Species = s__
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

# Hardcoded levels for UNITE — only k__Fungi exists in database
TAXO_FIXED <- list(
    Domain     = "Eukaryota",
    Supergroup = "Obazoa",
    Division   = NA   # filled from assignTaxonomy Phylum (p__)
)

# UNITE prefix stripping — applied to assignTaxonomy output
# e.g. "p__Ascomycota" -> "Ascomycota"
TAXO_STRIP_PREFIX <- TRUE

# --- Output filenames --------------------------------------------------------
OUT_FASTA          <- file.path(STEP4_DIR, paste0(MARKER, "_Fasta.fasta"))
OUT_MR             <- file.path(STEP4_DIR, paste0(MARKER, "_MR.csv"))
OUT_ASS            <- file.path(STEP5_DIR, paste0(MARKER, "_ASS.csv"))
OUT_SEQTAB_RDS     <- file.path(STEP4_DIR, "seqtab_nochim.rds")
OUT_VSEARCH_BLAST6 <- file.path(STEP5_DIR, paste0(MARKER, "_vsearch.blast6"))
