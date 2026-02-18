# Results Directory Guide

This directory stores generated outputs from the pipeline.

## Subdirectories

- `results/tables/` - CSV outputs from each analysis step
- `results/figures/` - PNG figure outputs

## Expected Outputs by Step

### Step 02 (QC)

- `tables/qc_lof_scores.csv`
- `tables/qc_shapiro_normality.csv`
- `figures/qc_lof_score_distribution.png`
- `figures/qc_shapiro_pvalue_distribution.png`

### Step 03 (Main statistical analysis + ANCOVA)

- `tables/stats_lipid_all.csv`
- `tables/stats_lipid_male.csv`
- `tables/stats_lipid_female.csv`
- `tables/stats_category_all.csv`
- `tables/stats_category_male.csv`
- `tables/stats_category_female.csv`
- `tables/ancova_sex_lipid.csv`
- `tables/ancova_sex_category.csv`

### Step 04 (No-AD sensitivity)

- `tables/sensitivity_noad_dataset.csv`
- `tables/sensitivity_noad_lipid_*.csv`
- `tables/sensitivity_noad_category_*.csv`
- `tables/sensitivity_noad_lipid_zscore_*.csv`

### Step 05 (Visualization)

- `figures/volcano_all.png`
- `figures/volcano_male.png`
- `figures/volcano_female.png`
- `figures/volcano_sex_interaction.png`
- `figures/category_effects_all.png`
- `figures/lipid_distribution_plots/*.png`

## Reproducibility Note

`results/` contents are generated artifacts and are intentionally gitignored by default. Regenerate them by rerunning the scripts in `scripts/README.md`.
