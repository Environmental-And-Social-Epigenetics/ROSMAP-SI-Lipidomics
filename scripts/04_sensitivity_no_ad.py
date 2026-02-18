"""Step 04: Sensitivity analysis excluding AD pathology cases (niareagansc > 2)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import (  # noqa: E402
    DATA_PROCESSED_DIR,
    FINAL_FORMATTED_FILENAME,
    TABLES_DIR,
    ensure_project_dirs,
)
from src.data_utils import get_lipid_columns  # noqa: E402
from src.stats_utils import (  # noqa: E402
    RegressionSpec,
    compute_category_means,
    filter_no_ad,
    run_per_lipid_regression,
    split_by_sex,
    zscore_columns,
)


def run_no_ad_sensitivity(
    input_csv: Path = DATA_PROCESSED_DIR / FINAL_FORMATTED_FILENAME,
    output_tables_dir: Path = TABLES_DIR,
    include_pmi: bool = False,
) -> dict[str, Path]:
    """
    Run no-AD sensitivity models.

    - Filter dataset to `niareagansc > 2`
    - Refit per-lipid models on all/male/female cohorts
    - Run category-mean models
    - Run z-score scaled lipid sensitivity models
    """
    ensure_project_dirs()
    output_tables_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_csv}. Run scripts/01_data_processing.py first."
        )

    df = pd.read_csv(input_csv)
    no_ad_df = filter_no_ad(df, reagan_column="niareagansc", threshold=2)
    no_ad_path = output_tables_dir / "sensitivity_noad_dataset.csv"
    no_ad_df.to_csv(no_ad_path, index=False)

    lipid_cols = get_lipid_columns(no_ad_df)
    if not lipid_cols:
        raise ValueError("No lipid columns detected for no-AD sensitivity models.")

    predictors = ["SI_avg", "niareagansc", "age_death"]
    if include_pmi and "pmi" in no_ad_df.columns:
        predictors.append("pmi")

    spec = RegressionSpec(
        predictors=tuple(predictors),
        primary_predictor="SI_avg",
        min_n=20,
    )

    output_paths: dict[str, Path] = {"noad_dataset": no_ad_path}
    cohorts = split_by_sex(no_ad_df)

    for cohort_name, cohort_df in cohorts.items():
        lipid_results = run_per_lipid_regression(cohort_df, lipid_columns=lipid_cols, spec=spec)
        lipid_out = output_tables_dir / f"sensitivity_noad_lipid_{cohort_name}.csv"
        lipid_results.to_csv(lipid_out, index=False)
        output_paths[f"lipid_{cohort_name}"] = lipid_out

        category_df = compute_category_means(cohort_df, lipid_columns=lipid_cols)
        category_cols = [c for c in category_df.columns if c.startswith("catmean_")]
        category_results = run_per_lipid_regression(
            category_df,
            lipid_columns=category_cols,
            spec=spec,
        )
        category_out = output_tables_dir / f"sensitivity_noad_category_{cohort_name}.csv"
        category_results.to_csv(category_out, index=False)
        output_paths[f"category_{cohort_name}"] = category_out

        z_df = zscore_columns(cohort_df, lipid_cols)
        z_results = run_per_lipid_regression(z_df, lipid_columns=lipid_cols, spec=spec)
        z_out = output_tables_dir / f"sensitivity_noad_lipid_zscore_{cohort_name}.csv"
        z_results.to_csv(z_out, index=False)
        output_paths[f"lipid_zscore_{cohort_name}"] = z_out

    print("No-AD sensitivity analysis complete.")
    for key, path in sorted(output_paths.items()):
        rows = len(pd.read_csv(path)) if path.exists() else 0
        print(f"{key}: {path} | rows={rows}")
    return output_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run no-AD lipidomics sensitivity analyses.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DATA_PROCESSED_DIR / FINAL_FORMATTED_FILENAME,
        help="Processed dataset to filter and model.",
    )
    parser.add_argument(
        "--output-tables-dir",
        type=Path,
        default=TABLES_DIR,
        help="Directory for sensitivity analysis outputs.",
    )
    parser.add_argument(
        "--include-pmi",
        action="store_true",
        help="Include pmi as an additional predictor when available.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_no_ad_sensitivity(
        input_csv=args.input_csv,
        output_tables_dir=args.output_tables_dir,
        include_pmi=args.include_pmi,
    )


if __name__ == "__main__":
    main()
