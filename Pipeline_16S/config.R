# =============================================================================
# File      : config.R
# Project   : MarkerFlow – generalised metabarcoding pipeline
# Marker    : 16S (Bacteria/Archaea)
# Author    : David Singer
# Date      : 2026-03-03
# Description: DADA2 and taxonomy parameters for 16S marker.
#              Sourced by 04_DADA2_ASV.R and 05_taxonomy.R.
#              PROJECT_ROOT is derived from this file's path in the R scripts.
# =============================================================================

# --- Marker ------------------------------------------------------------------
MARKER            <- "16S"
TAXONOMIC_GROUP   <- "bacteria"

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
DB_SILVA_VSEARCH  <- file.path(.project_root, "databases", "SILVA", "SILVA_138.2_SSURef_NR99_tax_silva.fasta")
DB_SILVA_DADA2    <- file.path(.project_root, "databases", "SILVA", "silva_nr99_v138.2_train_set.fa.gz")

# --- File naming -------------------------------------------------------------
FASTQ_STRIP_PATTERN_R <- "_R1\\.fastq\\.gz$"   # regex to strip R1 suffix → sample name

# --- DADA2 filterAndTrim parameters ------------------------------------------
# NOTE: Set TRUNCLEN based on read length distribution from MultiQC
# 16S amplicon is fixed-length — truncation is safe
# Adjust these values after checking trimmed read quality profiles
TRUNCLEN_FWD      <- 275
TRUNCLEN_REV      <- 275
MAX_EE_FWD        <- 1
MAX_EE_REV        <- 1
TRUNC_Q           <- 2
MAX_N             <- 0
COMPRESS          <- TRUE

# --- DADA2 denoising parameters ----------------------------------------------
POOLING           <- FALSE
CHIMERA_METHOD    <- "consensus"
THREADS           <- 32
THREADS_LEARN     <- 1            # learnErrors only — keep at 1 to avoid Linux segfault
MIN_OVERLAP       <- 12           # 16S V4-V5 ~411bp fixed amplicon — default overlap

# --- VSEARCH parameters ------------------------------------------------------
VSEARCH_PATH      <- "vsearch"    # full path if not in PATH
VSEARCH_ID        <- 0.1
VSEARCH_MAXACCEPTS <- 1
VSEARCH_MAXREJECTS <- 32
VSEARCH_THREADS   <- 32

# --- Taxonomy levels ---------------------------------------------------------
# SILVA 138.2 has 6 native levels mapped to 9-level standard:
#   Domain     = Domain    (Bacteria / Archaea / Eukaryota)
#   Supergroup = Phylum
#   Division   = NA        (not available in SILVA)
#   Subdivision= NA        (not available in SILVA)
#   Class      = Class
#   Order      = Order
#   Family     = Family
#   Genus      = Genus
#   Species    = NA        (not available in SILVA NR99)
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

# SILVA levels returned by assignTaxonomy — used for column mapping
# assignTaxonomy returns: Kingdom, Phylum, Class, Order, Family, Genus
# These map to our standard as: Domain, Supergroup, Class, Order, Family, Genus
# Division, Subdivision, Species are filled with NA
TAXO_SILVA_LEVELS <- c("Domain", "Supergroup", "Class", "Order", "Family", "Genus")
TAXO_SILVA_MAP    <- c(
    "Domain"      = "Domain",
    "Supergroup"  = "Supergroup",
    "Division"    = NA,
    "Subdivision" = NA,
    "Class"       = "Class",
    "Order"       = "Order",
    "Family"      = "Family",
    "Genus"       = "Genus",
    "Species"     = NA
)

# No hardcoded fixed levels for SILVA — multiple domains present
TAXO_FIXED        <- NULL
TAXO_STRIP_PREFIX <- FALSE

# --- Output filenames --------------------------------------------------------
OUT_FASTA          <- file.path(STEP4_DIR, paste0(MARKER, "_Fasta.fasta"))
OUT_MR             <- file.path(STEP4_DIR, paste0(MARKER, "_MR.csv"))
OUT_ASS            <- file.path(STEP5_DIR, paste0(MARKER, "_ASS.csv"))
OUT_SEQTAB_RDS     <- file.path(STEP4_DIR, "seqtab_nochim.rds")
OUT_VSEARCH_BLAST6 <- file.path(STEP5_DIR, paste0(MARKER, "_vsearch.blast6"))
