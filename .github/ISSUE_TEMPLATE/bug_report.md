---
name: Bug report
about: Report a problem with the pipeline
title: "[BUG] "
labels: bug
assignees: ''
---

**Describe the bug**
A clear and concise description of what went wrong.

**Command run**
```bash
# e.g. bash scripts/main.sh Pipeline_18S/config.sh --steps 3
```

**Marker / configuration**
- Marker: (18S / ITS / 16S / other)
- Relevant `config.sh` / `config.R` settings:

**Error message / logs**
```
Paste the error and the matching lines from results/<MARKER>/logs/
```

**Environment**
- OS: (e.g. Ubuntu 22.04 on WSL2)
- cutadapt: `cutadapt --version`
- vsearch: `vsearch --version`
- R / DADA2: `Rscript -e 'packageVersion("dada2")'`

**Additional context**
Anything else that might help (input data layout, sample naming, etc.).
