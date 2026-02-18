"""Step 05: Generate summary figures from statistical analysis outputs."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import (  # noqa: E402
    DATA_PROCESSED_DIR,
    FIGURES_DIR,
    FINAL_FORMATTED_FILENAME,
    TABLES_DIR,
    ensure_project_dirs,
)


def _ensure_stats_tables() -> None:
    """Ensure required stats tables exist before plotting."""
    required = [
        TABLES_DIR / "stats_lipid_all.csv",
        TABLES_DIR / "stats_lipid_male.csv",
        TABLES_DIR / "stats_lipid_female.csv",
        TABLES_DIR / "stats_category_all.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        missing_str = "\n".join(missing)
        raise FileNotFoundError(
            "Missing statistical analysis tables. Run scripts/03_statistical_analysis.py first.\n"
            f"{missing_str}"
        )


def make_volcano_plot(stats_df: pd.DataFrame, title: str, output_path: Path) -> None:
    """Create a volcano-style plot from per-lipid regression results."""
    if stats_df.empty:
        return

    plot_df = stats_df.copy()
    plot_df["neg_log10_p"] = -plot_df["p_value_primary"].clip(lower=1e-300).map(math.log10)
    plot_df["significant"] = plot_df["fdr_p_value"] < 0.05

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=plot_df,
        x="coef_primary",
        y="neg_log10_p",
        hue="significant",
        palette={True: "#d62728", False: "#1f77b4"},
        s=35,
        alpha=0.8,
        linewidth=0,
    )
    plt.axhline(-math.log10(0.05), linestyle="--", color="black", linewidth=1)
    plt.xlabel("SI_avg coefficient")
    plt.ylabel("-log10(p-value)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def make_category_barplot(category_df: pd.DataFrame, output_path: Path) -> None:
    """Plot category-level SI effect sizes (all cohort)."""
    if category_df.empty:
        return

    plot_df = category_df.copy()
    plot_df["category"] = plot_df["lipid"].str.replace("catmean_", "", regex=False)
    plot_df = plot_df.sort_values("coef_primary", ascending=False)

    plt.figure(figsize=(10, 6))
    sns.barplot(data=plot_df, x="category", y="coef_primary", color="#4c72b0")
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Lipid category")
    plt.ylabel("SI_avg coefficient")
    plt.title("Category-level SI association (all cohort)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def make_top_lipid_distribution_plots(
    data_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    output_dir: Path,
    top_n: int = 12,
) -> int:
    """Create scatter+trend plots for top-N lipids by p-value."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if stats_df.empty:
        return 0

    top = stats_df.sort_values("p_value_primary").head(top_n)
    created = 0

    for _, row in top.iterrows():
        lipid = row["lipid"]
        if lipid not in data_df.columns:
            continue

        plot_df = data_df[["SI_avg", lipid]].copy()
        plot_df = plot_df.dropna()
        if plot_df.empty:
            continue

        plt.figure(figsize=(6.5, 5))
        sns.regplot(
            data=plot_df,
            x="SI_avg",
            y=lipid,
            scatter_kws={"alpha": 0.65, "s": 20},
            line_kws={"color": "#d62728", "linewidth": 1.8},
            ci=None,
        )
        plt.title(f"{lipid}\ncoef={row['coef_primary']:.4f}, p={row['p_value_primary']:.2e}")
        plt.xlabel("SI_avg")
        plt.ylabel("Lipid level")
        plt.tight_layout()

        safe_name = lipid.replace("/", "_").replace(" ", "_")
        plt.savefig(output_dir / f"{safe_name}_vs_SI_avg.png", dpi=220)
        plt.close()
        created += 1

    return created


def run_visualization(
    processed_csv: Path = DATA_PROCESSED_DIR / FINAL_FORMATTED_FILENAME,
    tables_dir: Path = TABLES_DIR,
    figures_dir: Path = FIGURES_DIR,
    top_n_lipids: int = 12,
) -> dict[str, Path]:
    """Generate all summary visualizations."""
    ensure_project_dirs()
    _ensure_stats_tables()
    figures_dir.mkdir(parents=True, exist_ok=True)

    if not processed_csv.exists():
        raise FileNotFoundError(
            f"Processed dataset not found: {processed_csv}. Run scripts/01_data_processing.py first."
        )

    df = pd.read_csv(processed_csv)
    stats_all = pd.read_csv(tables_dir / "stats_lipid_all.csv")
    stats_male = pd.read_csv(tables_dir / "stats_lipid_male.csv")
    stats_female = pd.read_csv(tables_dir / "stats_lipid_female.csv")
    category_all = pd.read_csv(tables_dir / "stats_category_all.csv")

    outputs = {
        "volcano_all": figures_dir / "volcano_all.png",
        "volcano_male": figures_dir / "volcano_male.png",
        "volcano_female": figures_dir / "volcano_female.png",
        "category_barplot_all": figures_dir / "category_effects_all.png",
    }

    make_volcano_plot(stats_all, "Volcano plot: all cohort", outputs["volcano_all"])
    make_volcano_plot(stats_male, "Volcano plot: male cohort", outputs["volcano_male"])
    make_volcano_plot(stats_female, "Volcano plot: female cohort", outputs["volcano_female"])
    make_category_barplot(category_all, outputs["category_barplot_all"])

    dist_dir = figures_dir / "lipid_distribution_plots"
    count = make_top_lipid_distribution_plots(
        data_df=df,
        stats_df=stats_all,
        output_dir=dist_dir,
        top_n=top_n_lipids,
    )

    print("Visualization generation complete.")
    for key, path in outputs.items():
        print(f"{key}: {path}")
    print(f"top_lipid_distribution_plots: {dist_dir} | files_created={count}")
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate analysis summary visualizations.")
    parser.add_argument(
        "--processed-csv",
        type=Path,
        default=DATA_PROCESSED_DIR / FINAL_FORMATTED_FILENAME,
        help="Processed dataset used for scatter/distribution plots.",
    )
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=TABLES_DIR,
        help="Directory containing statistical result tables.",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=FIGURES_DIR,
        help="Directory for figure outputs.",
    )
    parser.add_argument(
        "--top-n-lipids",
        type=int,
        default=12,
        help="Number of top lipids (by p-value) to plot against SI_avg.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_visualization(
        processed_csv=args.processed_csv,
        tables_dir=args.tables_dir,
        figures_dir=args.figures_dir,
        top_n_lipids=args.top_n_lipids,
    )


if __name__ == "__main__":
    main()
