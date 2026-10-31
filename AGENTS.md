# Highlander package guide

## Package structure

- `R/` contains the package functions. `Highlander()` coordinates CMA-ES and
  LaplacesDemon; `Lowlander()` performs an initial Latin-hypercube search;
  `Lowlander_Z.R` contains helpers for reusing Lowlander results; and
  `Specs_help()` builds algorithm-specific LaplacesDemon settings.
- `man/` contains the package reference documentation, and `vignettes/` contains
  longer usage examples.
- `tests/testthat/` contains the package's unit tests. Testthat edition 3 is
  configured in `DESCRIPTION`.

## Tests and checks

- Run the unit tests from the repository root with
  `Rscript -e 'testthat::test_local()'`.
- Run the full package check with `R CMD build .` followed by
  `R CMD check --no-manual --as-cran Highlander_*.tar.gz`.
- Keep unit tests deterministic and fast. Prefer supplied Latin-hypercube
  matrices and small function inputs over random designs or full sampler runs.
- The `R-CMD-check` GitHub Actions workflow runs `R CMD check`, which includes
  the testthat suite. Add separate CI test steps only if they provide coverage
  that the package check does not.

## Change guidance

- Keep user-facing function behavior and argument defaults aligned with the
  reference documentation and README examples.
- When changing exported functions, update their reference documentation and
  add or adjust focused testthat coverage.
- Avoid tests that require a long MCMC/CMA optimization run when the same
  behavior can be checked with a small deterministic unit test.
