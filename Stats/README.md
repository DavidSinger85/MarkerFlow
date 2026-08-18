# MarkerFlow — Statistics stage

A small, standalone R workflow that turns the pipeline's per-marker outputs into filtered
community matrices and a set of dataset-quality diagnostic figures. It is **marker-agnostic** and
**metadata-free** — it needs only the `MR` (read matrix) and `ASS` (taxonomy) files produced by the
bioinformatics pipeline.

```
00_import.R          → import MR + ASS per marker  → <PROJECT>.RData (object DS)
01_filtering.R       → taxonomic / rare-ASV / quality / sample filtering,
                       averaged rarefaction, read-depth residual correction
02_dataset_quality.R → diagnostic figures + summary tables
```

## Prerequisites

1. **Run the bioinformatics pipeline first.** This stage reads
   `../results/<MARKER>/04_ASV/<MARKER>_MR.csv` and
   `../results/<MARKER>/05_taxonomy/<MARKER>_ASS.csv`. Those live under the gitignored `results/`
   directory, so they are **not** in a fresh clone — generate them with
   `bash scripts/main.sh Pipeline_<MARKER>/config.sh --steps 1,2,3,4,5`.
2. **R packages:** `vegan`, `ggplot2`, `patchwork`, `scales`. Install once, e.g. with conda
   (`conda install -c conda-forge r-vegan r-patchwork r-ggplot2 r-scales`) or in R
   (`install.packages(c("vegan","ggplot2","patchwork","scales"))`).

## Running it

Open `Stats.Rproj` in RStudio (this sets the working directory to `Stats/`), then run the three
scripts in order:

```r
source("00_import.R")
source("01_filtering.R")
source("02_dataset_quality.R")
```

Or from a shell:

```bash
cd Stats
Rscript 00_import.R && Rscript 01_filtering.R && Rscript 02_dataset_quality.R
```

## Configuration

Each script has a **Section 1 — Parameters** block, the single place to edit. Key settings:

- `PROJECT_NAME` — names the `.RData` object file and all outputs (`00`).
- `MARKERS`, `MARKER_LABELS`, `MARKER_COLORS` — which markers to process and how to label them.
- `RESULTS_DIR` — where the pipeline outputs live (defaults to `../results`).
- `EXCLUDE_TAXA` — optional per-marker taxonomic exclusion (empty by default; a PR2/18S "protists
  only" example is shown commented in `01`).
- `EXCLUDE_SAMPLES`, `FILTER_MIN_*` — sample/ASV filtering thresholds.
- `RARE_THRESHOLDS` — per-marker rarefaction depth. **The defaults are low, sized for the small
  bundled example data** — raise them for a real dataset (guided by the read-depth figures produced
  by `01`/`02`).

## Outputs

- `Stats/<PROJECT>.RData` — the `DS` object. Per marker: `MR`, `ASS`, `MR_clean`, `ASS_clean`,
  `MR_rare`, `ASS_rare`, `MR_cor`, `PARAMS`.
- `Stats/outputs/figures/` — diagnostic figures (PDF + JPEG): read-depth distribution, reads vs
  observed richness across pipeline stages, rarefaction curves, ASV length distribution, % identity
  distribution, filtering journey, and the decision-support figures from `01`.
- `Stats/outputs/tables/` — filtering-journey and read-depth summary CSVs.

`outputs/` and `<PROJECT>.RData` are gitignored (regenerated on every run).
