# Contributing to MarkerFlow

Thanks for your interest in improving MarkerFlow! Contributions of all kinds are
welcome: bug reports, documentation fixes, new marker configurations, and code.

## Reporting bugs and requesting features

Please open an issue using the templates provided (bug report / feature request).
For bugs, include:

- The command you ran (e.g. `bash scripts/main.sh Pipeline_18S/config.sh --steps 3`)
- The marker and the relevant part of your `config.sh` / `config.R`
- The error message and, if useful, the matching lines from
  `results/<MARKER>/logs/`
- Tool versions (`cutadapt --version`, `Rscript -e 'packageVersion("dada2")'`, etc.)

## Proposing a new marker

MarkerFlow is marker-agnostic — most new markers need only a new configuration,
not new code:

1. Copy the templates:
   ```bash
   cp config/config_template.sh Pipeline_<MARKER>/config.sh
   cp config/config_template.R  Pipeline_<MARKER>/config.R
   ```
2. Fill in the primers, run directories, database paths, and DADA2 parameters.
3. Add a small example dataset under `RAW_dataset/<MARKER>/` if you can share one.
4. Open a pull request describing the marker and the reference database used.

## Submitting changes

1. Fork the repository and create a branch off `main`.
2. Keep changes focused; match the existing style of the surrounding code.
3. Before opening a PR, run the same checks the CI runs:
   ```bash
   shellcheck scripts/*.sh Pipeline_*/config.sh
   Rscript -e 'invisible(lapply(list.files(c("scripts","Pipeline_18S"), pattern="[.]R$", full.names=TRUE), parse))'
   ```
4. Update `README.md` and `CHANGELOG.md` if your change affects users.
5. Open a pull request using the template and link any related issue.

## Coding conventions

- **Bash**: `#!/usr/bin/env bash`, `set -euo pipefail`, functions from `scripts/utils.sh`
  for logging and read tracking. No hardcoded paths — everything derives from
  `PROJECT_ROOT` / the config file location.
- **R**: parameters live in `config.R`, not in the analysis scripts. Keep the
  9-level taxonomy standard (`Domain … Species`) intact.

## Code of conduct

By participating you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).
