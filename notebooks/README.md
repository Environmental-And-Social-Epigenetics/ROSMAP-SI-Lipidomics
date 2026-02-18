# Notebook Workflow Guide

This directory contains the interactive, inline-code walkthrough of the full lipidomics analysis.

## Run Order

Run notebooks in this sequence:

1. `01_data_processing.ipynb`
2. `02_quality_control.ipynb`
3. `03_statistical_analysis.ipynb`
4. `04_sensitivity_no_ad.ipynb`
5. `05_visualization.ipynb`

## Important Run Notes

- Start Jupyter from the **repository root** so imports like `from config import ...` resolve correctly.
- Each notebook is self-contained and does not rely on `subprocess.run(...)`.
- Notebooks write outputs into `data/processed/` and `results/`.

## Notebook Summary

| Notebook | Purpose | Key inputs | Key outputs |
|---|---|---|---|
| `01_data_processing.ipynb` | Build analysis-ready sample-by-lipid datasets | raw CSVs in `data/raw/` | `data/processed/Normalized_Formatted_Lipidomics.csv`, `data/processed/Final_Formatted_Lipidomics.csv` |
| `02_quality_control.ipynb` | Outlier and normality diagnostics | `Final_Formatted_Lipidomics.csv` | `results/tables/qc_*.csv`, QC histograms |
| `03_statistical_analysis.ipynb` | Core OLS + ANCOVA interaction modeling | `Final_Formatted_Lipidomics.csv` | `results/tables/stats_*`, `ancova_sex_*` |
| `04_sensitivity_no_ad.ipynb` | Repeat analysis in no-AD cohort (`niareagansc > 2`) | `Final_Formatted_Lipidomics.csv` | `results/tables/sensitivity_noad_*` |
| `05_visualization.ipynb` | Generate summary figures including ANCOVA volcano | stats/sensitivity tables + processed dataset | `results/figures/*.png` |

## Recommended Usage Pattern

- Use scripts for reproducible batch runs.
- Use notebooks for exploration, interpretation, and figure review.
