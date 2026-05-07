# R2easyviz (development version)

## Bug fixes and maintenance

* Fixed `Author@R` typo in `DESCRIPTION` (now correctly `Authors@R`).
* Declared previously-undeclared imports (`pheatmap`, `viridis`, `rlang`,
  `tidyr`) so installation succeeds for users without those packages
  already loaded.
* Replaced deprecated `slot = "data"` argument with `layer = "data"` in
  `r2correlate_samples()` for compatibility with Seurat v5.
* Corrected `@return` documentation in `r2correlate_samples()` (the
  matrix branch was previously labelled `return_heatmap = TRUE`).
* Added missing `@export` tag to `r2prop_plot()`.

## Internal

* Switched broad `@import` directives to specific `@importFrom` calls
  and namespace-qualified all external function calls.
* Replaced wildcard `exportPattern` in `NAMESPACE` with explicit
  roxygen2-managed exports.
* Removed scaffolding `R/hello.R` and dev artifacts (`save.RData`,
  `testing.Rmd`, `filtered_gene_bc_matrices/`) from the repo.
* Added `testthat` test scaffolding under `tests/testthat/`.
* Added GitHub Actions `R-CMD-check` workflow.

# R2easyviz 0.1.0

* Initial development release.
