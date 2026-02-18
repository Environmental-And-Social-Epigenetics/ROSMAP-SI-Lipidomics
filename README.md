# ROSMAP-SI-Lipidomics

Rebuilt, fully documented analysis pipeline for ROSMAP social isolation and brain lipidomics associations.

## Project Goal

This repository tests whether social isolation (primarily `SI_avg`) is associated with lipid abundance in ROSMAP brain tissue, while adjusting for key covariates (for example `niareagansc` and `age_death`) and stratifying by sex.

The original exploratory analysis lived in:

- `/om2/user/mabdel03/files/Isolation/Lipidomics_Work/`

This repository reorganizes that work into a reproducible, collaborator-friendly structure with clear script and notebook entrypoints.

## What This Rebuild Preserves

- Data-processing logic from `SI_Lipidomics_Processing.ipynb`
- QC logic from `Lipidomics_Code_Review.ipynb`
- Full-cohort/sex-stratified modeling ideas from `R_stats.ipynb`, `Isolation_Lipidomics_Stats.ipynb`, and `Rebecca_Replication_Analysis.ipynb`
- No-AD sensitivity analyses from `Lipidomics_NoAD.ipynb` and `Lipidomics_NoAD_March2025Update.ipynb`
- Figure-generation concepts from `Lipidomics_Final_Figures.ipynb`

## Repository Structure

```text
ROSMAP-SI-Lipidomics/
├── config.py
├── requirements.txt
├── data/
│   ├── README.md
│   ├── raw/                      # local copies of source CSV files (gitignored)
│   └── processed/                # pipeline outputs (gitignored)
├── src/
│   ├── data_utils.py             # data loading/prep helpers
│   └── stats_utils.py            # QC + model helper functions
├── scripts/
│   ├── 01_data_processing.py
│   ├── 02_quality_control.py
│   ├── 03_statistical_analysis.py
│   ├── 04_sensitivity_no_ad.py
│   └── 05_visualization.py
├── notebooks/
│   ├── 01_data_processing.ipynb
│   ├── 02_quality_control.ipynb
│   ├── 03_statistical_analysis.ipynb
│   ├── 04_sensitivity_no_ad.ipynb
│   └── 05_visualization.ipynb
└── results/
    ├── figures/                  # generated plots (gitignored)
    └── tables/                   # generated tables (gitignored)
```

## Data Handling

- Source CSV files are copied locally into `data/raw/` for reproducibility.
- Original files in `/om2/user/mabdel03/files/Isolation/Lipidomics_Work/` are not edited.
- CSVs are intentionally gitignored to keep the repository lightweight.
- Full file inventory and provenance are documented in `data/README.md`.

## Setup

From the repository root:

```bash
python -m pip install -r requirements.txt
```

## End-to-End Run Order

Run scripts in order:

```bash
python scripts/01_data_processing.py
python scripts/02_quality_control.py
python scripts/03_statistical_analysis.py
python scripts/04_sensitivity_no_ad.py
python scripts/05_visualization.py
```

## Pipeline Steps

### Step 01 - Data Processing

Script: `scripts/01_data_processing.py`  
Notebook: `notebooks/01_data_processing.ipynb`

- Loads raw lipidomics matrix and metadata
- Maps `projid -> individualID`
- Reshapes lipidomics matrix to sample-by-lipid format
- Applies internal standard normalization (`PC(18:1D7_15:0)` positive/negative references)
- Attaches covariates (`SI_avg`, `niareagansc`, `msex`, `age_death`, `educ`, medication vars, `pmi`)
- Writes:
  - `data/processed/Normalized_Formatted_Lipidomics.csv`
  - `data/processed/Final_Formatted_Lipidomics.csv`

### Step 02 - Quality Control

Script: `scripts/02_quality_control.py`  
Notebook: `notebooks/02_quality_control.ipynb`

- Local Outlier Factor (LOF) sample-level outlier scan
- Shapiro-Wilk normality test per lipid
- Benjamini-Hochberg FDR adjustment
- Writes:
  - `results/tables/qc_lof_scores.csv`
  - `results/tables/qc_shapiro_normality.csv`
  - `results/figures/qc_lof_score_distribution.png`
  - `results/figures/qc_shapiro_pvalue_distribution.png`

### Step 03 - Statistical Analysis (Full Cohort)

Script: `scripts/03_statistical_analysis.py`  
Notebook: `notebooks/03_statistical_analysis.ipynb`

- Per-lipid OLS model:
  - `lipid ~ SI_avg + niareagansc + age_death`
- Cohorts:
  - all
  - male (`msex == 1`)
  - female (`msex == 0`)
- Category-mean models (`catmean_*`)
- Writes:
  - `results/tables/stats_lipid_all.csv`
  - `results/tables/stats_lipid_male.csv`
  - `results/tables/stats_lipid_female.csv`
  - `results/tables/stats_category_all.csv`
  - `results/tables/stats_category_male.csv`
  - `results/tables/stats_category_female.csv`

### Step 04 - Sensitivity Analysis (No AD)

Script: `scripts/04_sensitivity_no_ad.py`  
Notebook: `notebooks/04_sensitivity_no_ad.ipynb`

- Filters to no/low AD pathology (`niareagansc > 2`)
- Re-runs per-lipid and category models for all/male/female cohorts
- Adds z-score scaled lipid sensitivity models
- Writes:
  - `results/tables/sensitivity_noad_dataset.csv`
  - `results/tables/sensitivity_noad_lipid_*.csv`
  - `results/tables/sensitivity_noad_category_*.csv`
  - `results/tables/sensitivity_noad_lipid_zscore_*.csv`

### Step 05 - Visualization

Script: `scripts/05_visualization.py`  
Notebook: `notebooks/05_visualization.ipynb`

- Creates:
  - Volcano plots (`all`, `male`, `female`)
  - Category-level effect-size barplot
  - Top lipid regression/distribution plots vs `SI_avg`
- Writes:
  - `results/figures/volcano_all.png`
  - `results/figures/volcano_male.png`
  - `results/figures/volcano_female.png`
  - `results/figures/category_effects_all.png`
  - `results/figures/lipid_distribution_plots/*.png`

## Notes for Collaborators

- Start with scripts for reproducible command-line runs.
- Use notebooks for interactive exploration and interpretation.
- Path constants are centralized in `config.py`.
- Shared preprocessing/model code is centralized under `src/`.
