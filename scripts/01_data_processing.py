"""Step 01: Build analysis-ready lipidomics datasets from raw input tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config import (
    DATA_PROCESSED_DIR,
    FINAL_FORMATTED_FILENAME,
    NORMALIZED_FORMATTED_FILENAME,
    ensure_project_dirs,
)
from src.data_utils import (
    add_individual_ids_to_longitudinal,
    attach_covariates,
    load_raw_input_tables,
    merge_lipid_id_map,
    normalize_by_internal_standards,
    reshape_lipidomics_matrix,
    add_individual_id_column,
    drop_optional_columns,
)


def build_outputs(
    output_dir: Path = DATA_PROCESSED_DIR,
    require_si_for_final: bool = True,
) -> tuple[Path, Path]:
    """
    Generate normalized and final formatted datasets for downstream analyses.

    Returns:
        Tuple of (normalized_output_path, final_output_path).
    """
    ensure_project_dirs()
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = load_raw_input_tables()

    longitudinal = add_individual_ids_to_longitudinal(
        tables["metadata_longitudinal"],
        tables["metadata_clinical"],
    )

    sample_df = reshape_lipidomics_matrix(tables["lipidomics_raw"])
    sample_df = add_individual_id_column(sample_df)
    sample_df = normalize_by_internal_standards(sample_df)
    sample_df = merge_lipid_id_map(sample_df, tables["lipid_individual_map"])
    sample_df = attach_covariates(sample_df, longitudinal)
    sample_df = drop_optional_columns(sample_df, columns=["Net_SI", "Sex"])

    normalized_path = output_dir / NORMALIZED_FORMATTED_FILENAME
    sample_df.to_csv(normalized_path, index=False)

    final_df = sample_df.copy()
    final_df = final_df.dropna(subset=["individualID"])
    if require_si_for_final and "SI_avg" in final_df.columns:
        final_df = final_df.dropna(subset=["SI_avg"])
    final_df = final_df.drop_duplicates(subset=["label"], keep="first")

    final_path = output_dir / FINAL_FORMATTED_FILENAME
    final_df.to_csv(final_path, index=False)

    print("Data processing complete.")
    print(f"Normalized dataset: {normalized_path} | shape={sample_df.shape}")
    print(f"Final dataset:      {final_path} | shape={final_df.shape}")
    return normalized_path, final_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create normalized/final formatted datasets for ROSMAP lipidomics analysis.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_PROCESSED_DIR,
        help="Destination directory for processed CSV outputs.",
    )
    parser.add_argument(
        "--allow-missing-si-final",
        action="store_true",
        help="If provided, keep rows with missing SI_avg in Final_Formatted output.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_outputs(
        output_dir=args.output_dir,
        require_si_for_final=not args.allow_missing_si_final,
    )


if __name__ == "__main__":
    main()
