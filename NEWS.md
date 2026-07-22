# DPM 0.1.0

## Initial public release

* Added the complete DPM workflow for polygenic risk decomposition using a heteroscedastic Dirichlet process mixture.
* Added single-cell expression aggregation with `preprocess_scrna()`.
* Added GWAS summary-statistics preprocessing with `preprocess_gwas()`.
* Added gene annotation readers and SNP-to-gene mapping through `read_gene_annotation()` and `match_snps_to_genes()`.
* Added LD-corrected gene-level response construction through `compute_response()`.
* Added SNP-name and chromosome-position matching modes.
* Added reciprocal completion and symmetrization of sparsely stored LD pairs.
* Added trace-preserving PSD repair through eigenvalue clipping.
* Added chromosome-level response output and `combine_response_by_chr()`.
* Added gene2vec and k-means initialization interfaces.
* Added prior recommendation and model-input construction through `recommend_priors()` and `prepare_input()`.
* Added MCMC fitting through `run_dpm()` with validation, progress reporting, retained-draw diagnostics, and optional scale storage.
* Added posterior-similarity-based partitioning through `get_partition()`.
* Added combined FGLS and Bayesian association inference through `test_significance()`.
* Added protection for underdetermined frequentist module fits.
* Added gene and cellular-label prioritization through `gene_score()` and `cell_type_score()`.
* Added a 200-gene installed example dataset for regression tests, integration tests, and documentation.
* Added three package vignettes covering the core workflow, response computation, and input formats.
* Added automated `testthat` coverage for example-data integrity, response reproducibility, and the end-to-end workflow.
* Verified package installation, examples, compiled code, tests, and vignette rebuilding through `R CMD check`.
