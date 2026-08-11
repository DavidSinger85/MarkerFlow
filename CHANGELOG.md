# Changelog

All notable changes to MetaBarFlow are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-08-10

First public release.

### Added
- Five-step metabarcoding pipeline: merge → quality check → primer trimming →
  DADA2 ASV inference → dual (VSEARCH + DADA2) taxonomy assignment.
- Marker-agnostic design driven entirely by per-marker configuration files
  (`Pipeline_<MARKER>/config.sh` + `config.R`).
- Ready-to-use configurations for three markers: 18S (PR2), ITS (UNITE), 16S (SILVA).
- Documented configuration templates in `config/` for adding new markers.
- Per-sample read tracking across every step (`results/<MARKER>/logs/read_tracking.tsv`).
- Small real example dataset for all three markers under `RAW_dataset/` so a fresh
  install can be verified end to end. FASTQ files are real reads subsampled to ~3,000 read
  pairs each (`test/subsample_raw.sh`) to keep the repository lightweight.
- Statistics stage under `Stats/` (`00_import` → `01_filtering` →
  `02_dataset_quality`): a standalone, marker-agnostic, metadata-free R workflow that imports the
  pipeline's `MR`/`ASS` outputs, applies taxonomic/rare/quality/sample filtering, averaged
  rarefaction and read-depth residual correction, and produces dataset-quality diagnostic figures
  and tables.
- `README.md`, `INSTALL.md`, and a test guide under `test/`.
- Standard open-source project files: LICENSE (MIT), CONTRIBUTING, CODE_OF_CONDUCT,
  CITATION, issue/PR templates, and a lint CI workflow.

### Changed
- Pipeline output moved from `pipeline_<MARKER>/` to `results/<MARKER>/` so that
  configuration directories (`Pipeline_<MARKER>/`) and output directories no longer
  collide on case-insensitive filesystems (Windows / WSL).
- Renamed the marker configuration directories from `Piepline_*` to `Pipeline_*`.

### Fixed
- DADA2 read-tracking columns (`step04_merged`, `step04_nochim`).
- `learnErrors` now uses a dedicated `THREADS_LEARN` (default 1) to avoid a Linux segfault.
- RDS checkpoints added after `learnErrors` and `dada` for crash recovery.
- Sample-name derivation is now config-driven, not hardcoded to a project prefix.
- `main.sh` usage examples now reference the correct configuration paths.

[1.0.0]: https://github.com/DavidSinger85/MetaBarFlow/releases/tag/v1.0.0
