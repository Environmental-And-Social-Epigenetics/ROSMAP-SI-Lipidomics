# Script Workflow Guide

This directory contains the command-line implementation of the full analysis pipeline.

## Sequential Run Order

```bash
python scripts/01_data_processing.py
python scripts/02_quality_control.py
python scripts/03_statistical_analysis.py
python scripts/04_sensitivity_no_ad.py
python scripts/05_visualization.py
```

## One-liner Full Run

```bash
python scripts/01_data_processing.py && \
python scripts/02_quality_control.py && \
python scripts/03_statistical_analysis.py && \
python scripts/04_sensitivity_no_ad.py && \
python scripts/05_visualization.py
```

## Script Summary

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `01_data_processing.py` | Build processed analysis datasets from raw tables | `data/raw/*` | `data/processed/Normalized_Formatted_Lipidomics.csv`, `data/processed/Final_Formatted_Lipidomics.csv` |
| `02_quality_control.py` | LOF outlier + Shapiro normality QC | `data/processed/Final_Formatted_Lipidomics.csv` | `results/tables/qc_*.csv`, QC plots |
| `03_statistical_analysis.py` | Main OLS models + ANCOVA sex-interaction models | `data/processed/Final_Formatted_Lipidomics.csv` | `results/tables/stats_*`, `results/tables/ancova_sex_*` |
| `04_sensitivity_no_ad.py` | No-AD and scaling sensitivity analyses | `data/processed/Final_Formatted_Lipidomics.csv` | `results/tables/sensitivity_noad_*` |
| `05_visualization.py` | Summary figures from model outputs | `results/tables/*`, processed data | `results/figures/*.png` |

## Common CLI Flags

- `scripts/03_statistical_analysis.py --include-pmi`
  - Adds `pmi` as an additional covariate when available.
- `scripts/04_sensitivity_no_ad.py --include-pmi`
  - Adds `pmi` in sensitivity models.
- `scripts/05_visualization.py --top-n-lipids 20`
  - Controls number of top lipids plotted in distribution/scatter outputs.
