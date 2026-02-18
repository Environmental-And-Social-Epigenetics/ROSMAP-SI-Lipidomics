"""Step 02: Run quality-control diagnostics for the processed lipidomics dataset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

from config import (  # noqa: E402
    DATA_PROCESSED_DIR,
    FIGURES_DIR,
    FINAL_FORMATTED_FILENAME,
    TABLES_DIR,
    ensure_project_dirs,
)
from src.data_utils import get_lipid_columns  # noqa: E402
from src.stats_utils import run_lof_outlier_detection, run_shapiro_normality_tests  # noqa: E402


def run_quality_control(
    input_csv: Path = DATA_PROCESSED_DIR / FINAL_FORMATTED_FILENAME,
    output_tables_dir: Path = TABLES_DIR,
    output_figures_dir: Path = FIGURES_DIR,
) -> tuple[Path, Path]:
    """Run LOF and normality diagnostics and save result tables/figures."""
    ensure_project_dirs()
    output_tables_dir.mkdir(parents=True, exist_ok=True)
    output_figures_dir.mkdir(parents=True, exist_ok=True)

    if not input_csv.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_csv}. Run scripts/01_data_processing.py first."
        )

    df = pd.read_csv(input_csv)
    lipid_cols = get_lipid_columns(df)
    if not lipid_cols:
        raise ValueError("No lipid feature columns were detected in the input dataset.")

    lof_df = run_lof_outlier_detection(df, lipid_columns=lipid_cols, n_neighbors=100)
    normality_df = run_shapiro_normality_tests(df, lipid_columns=lipid_cols)

    lof_path = output_tables_dir / "qc_lof_scores.csv"
    normality_path = output_tables_dir / "qc_shapiro_normality.csv"
    lof_df.to_csv(lof_path, index=False)
    normality_df.to_csv(normality_path, index=False)

    # Diagnostic plot: LOF score distribution
    if plt is not None and not lof_df.empty:
        plt.figure(figsize=(8, 5))
        plt.hist(lof_df["lof_score"], bins=30, edgecolor="black")
        plt.title("LOF Score Distribution")
        plt.xlabel("Negative Outlier Factor")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(output_figures_dir / "qc_lof_score_distribution.png", dpi=200)
        plt.close()

    # Diagnostic plot: normality p-value distribution
    if plt is not None and not normality_df.empty:
        plt.figure(figsize=(8, 5))
        plt.hist(normality_df["p_value"], bins=30, edgecolor="black")
        plt.title("Shapiro-Wilk P-value Distribution Across Lipids")
        plt.xlabel("P-value")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(output_figures_dir / "qc_shapiro_pvalue_distribution.png", dpi=200)
        plt.close()

    print("Quality-control analysis complete.")
    print(f"LOF table: {lof_path} | rows={len(lof_df)}")
    print(f"Normality table: {normality_path} | rows={len(normality_df)}")
    return lof_path, normality_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run QC diagnostics for lipidomics data.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DATA_PROCESSED_DIR / FINAL_FORMATTED_FILENAME,
        help="Processed dataset to evaluate.",
    )
    parser.add_argument(
        "--output-tables-dir",
        type=Path,
        default=TABLES_DIR,
        help="Directory for QC output tables.",
    )
    parser.add_argument(
        "--output-figures-dir",
        type=Path,
        default=FIGURES_DIR,
        help="Directory for QC diagnostic plots.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_quality_control(
        input_csv=args.input_csv,
        output_tables_dir=args.output_tables_dir,
        output_figures_dir=args.output_figures_dir,
    )


if __name__ == "__main__":
    main()
