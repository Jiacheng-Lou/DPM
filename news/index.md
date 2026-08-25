# Changelog

## DPM 0.1.0

### Initial public release

- Added the complete DPM workflow for polygenic risk decomposition using
  a heteroscedastic Dirichlet process mixture.
- Added single-cell expression aggregation with
  [`preprocess_scrna()`](https://jiacheng-lou.github.io/DPM/reference/preprocess_scrna.md).
- Added GWAS summary-statistics preprocessing with
  [`preprocess_gwas()`](https://jiacheng-lou.github.io/DPM/reference/preprocess_gwas.md).
- Added gene annotation readers and SNP-to-gene mapping through
  [`read_gene_annotation()`](https://jiacheng-lou.github.io/DPM/reference/read_gene_annotation.md)
  and
  [`match_snps_to_genes()`](https://jiacheng-lou.github.io/DPM/reference/match_snps_to_genes.md).
- Added LD-corrected gene-level response construction through
  [`compute_response()`](https://jiacheng-lou.github.io/DPM/reference/compute_response.md).
- Added SNP-name and chromosome-position matching modes.
- Added reciprocal completion and symmetrization of sparsely stored LD
  pairs.
- Added trace-preserving PSD repair through eigenvalue clipping.
- Added chromosome-level response output and
  [`combine_response_by_chr()`](https://jiacheng-lou.github.io/DPM/reference/combine_response_by_chr.md).
- Added gene2vec and k-means initialization interfaces.
- Added prior recommendation and model-input construction through
  [`recommend_priors()`](https://jiacheng-lou.github.io/DPM/reference/recommend_priors.md)
  and
  [`prepare_input()`](https://jiacheng-lou.github.io/DPM/reference/prepare_input.md).
- Added MCMC fitting through
  [`run_dpm()`](https://jiacheng-lou.github.io/DPM/reference/run_dpm.md)
  with validation, progress reporting, retained-draw diagnostics, and
  optional scale storage.
- Added posterior-similarity-based partitioning through
  [`get_partition()`](https://jiacheng-lou.github.io/DPM/reference/get_partition.md).
- Added combined FGLS and Bayesian association inference through
  [`test_significance()`](https://jiacheng-lou.github.io/DPM/reference/test_significance.md).
- Added protection for underdetermined frequentist module fits.
- Added gene and cellular-label prioritization through
  [`gene_score()`](https://jiacheng-lou.github.io/DPM/reference/gene_score.md)
  and
  [`cell_type_score()`](https://jiacheng-lou.github.io/DPM/reference/cell_type_score.md).
- Added a 200-gene installed example dataset for regression tests,
  integration tests, and documentation.
- Added three package vignettes covering the core workflow, response
  computation, and input formats.
- Added automated `testthat` coverage for example-data integrity,
  response reproducibility, and the end-to-end workflow.
- Verified package installation, examples, compiled code, tests, and
  vignette rebuilding through `R CMD check`.
