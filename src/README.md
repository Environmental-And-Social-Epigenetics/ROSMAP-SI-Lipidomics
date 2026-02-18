# Source Utilities Guide

This directory contains shared utility code used by both scripts and notebooks.

## Modules

### `data_utils.py`

Data ingestion and preprocessing helpers.

Key functions:

- `load_raw_input_tables()` - load all expected raw input CSVs.
- `add_individual_ids_to_longitudinal()` - map `projid` to `individualID`.
- `reshape_lipidomics_matrix()` - convert lipid-row matrix to sample-row matrix.
- `add_individual_id_column()` - parse `individualID` from sample labels.
- `normalize_by_internal_standards()` - normalize positive/negative mode lipids by internal standards.
- `merge_lipid_id_map()` - merge specimen-to-individual mapping.
- `attach_covariates()` - add analysis covariates (SI, pathology, demographics, meds, pmi).
- `get_lipid_columns()` - identify lipid feature columns.
- `prepare_analysis_dataset()` - end-to-end preprocessing helper.

### `stats_utils.py`

Statistical and QC helpers.

Key functions:

- `run_lof_outlier_detection()` - sample-level LOF outlier diagnostics.
- `run_shapiro_normality_tests()` - per-lipid normality tests.
- `add_fdr_column()` - Benjamini-Hochberg correction utility.
- `run_per_lipid_regression()` - batch OLS across lipids.
- `run_ancova_sex_interaction()` - pooled ANCOVA interaction model (`SI_avg * msex`).
- `compute_category_means()` - category-level mean features (`catmean_*`).
- `split_by_sex()` - all/male/female cohort split.
- `filter_no_ad()` - no-AD subset (`niareagansc > 2`).
- `zscore_columns()` - z-score scaling utility.

## Import Pattern

Scripts and notebooks import from `src` via the repository root. The scripts ensure this by appending project root to `sys.path` at runtime.
