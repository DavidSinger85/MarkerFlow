## Summary

Briefly describe what this pull request changes and why.

Closes # (issue number, if applicable)

## Type of change

- [ ] Bug fix
- [ ] New marker configuration
- [ ] New feature
- [ ] Documentation
- [ ] Other:

## Checklist

- [ ] I ran `shellcheck scripts/*.sh Pipeline_*/config.sh`
- [ ] R scripts still parse (`Rscript -e 'lapply(list.files("scripts", "[.]R$", full.names=TRUE), parse)'`)
- [ ] No hardcoded paths, secrets, or large data files added
- [ ] Updated `README.md` / `CHANGELOG.md` if user-facing behavior changed
- [ ] Tested on the example dataset (if applicable)

## Notes for reviewers

Anything reviewers should pay particular attention to.
