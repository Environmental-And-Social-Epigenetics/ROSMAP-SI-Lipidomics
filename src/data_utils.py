"""Data loading and preprocessing helpers for the lipidomics analysis pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from config import (
    DATA_RAW_DIR,
    LIPID_ID_MAP_FILENAME,
    META_CLINICAL_FILENAME,
    META_LONGITUDINAL_FILENAME,
    ORIGINAL_FINAL_FORMATTED_FILENAME,
    ORIGINAL_NORMALIZED_FORMATTED_FILENAME,
    RAW_LIPIDOMICS_FILENAME,
)

# Lipid categories used in the original notebooks.
LIPID_CATEGORY_SUFFIXES = {
    "AcCa",
    "Cer",
    "ChE",
    "DAG",
    "GSL",
    "LPC",
    "LPE",
    "MAG",
    "PC",
    "PE",
    "PG",
    "PI",
    "PS",
    "SM",
    "TG",
}

DEFAULT_COVARIATE_COLUMNS = [
    "social_isolation_avg",
    "social_isolation_lv",
    "niareagansc",
    "msex",
    "age_death",
    "educ",
    "lipid_lowering_nonstatin_rx_ever",
    "statin_rx_bl",
    "statin_rx_ever",
    "pmi",
]

RENAME_COVARIATES = {
    "social_isolation_avg": "SI_avg",
    "social_isolation_lv": "SI_lv",
}


def load_raw_input_tables(data_raw_dir: Path = DATA_RAW_DIR) -> dict[str, pd.DataFrame]:
    """Load all source input tables used by the rebuilt analysis pipeline."""
    paths = {
        "lipidomics_raw": data_raw_dir / "lipidomics_data" / RAW_LIPIDOMICS_FILENAME,
        "metadata_longitudinal": data_raw_dir / "metadata" / META_LONGITUDINAL_FILENAME,
        "metadata_clinical": data_raw_dir / "metadata" / META_CLINICAL_FILENAME,
        "lipid_individual_map": data_raw_dir / "metadata" / LIPID_ID_MAP_FILENAME,
        "original_final_formatted": data_raw_dir
        / "original_processed_csvs"
        / ORIGINAL_FINAL_FORMATTED_FILENAME,
        "original_normalized_formatted": data_raw_dir
        / "original_processed_csvs"
        / ORIGINAL_NORMALIZED_FORMATTED_FILENAME,
    }
    return {name: pd.read_csv(path) for name, path in paths.items()}


def add_individual_ids_to_longitudinal(
    longitudinal_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Attach ``individualID`` to the longitudinal table via ``projid``.

    This mirrors the mapping logic in the original preprocessing notebook.
    """
    mapping = (
        clinical_df[["projid", "individualID"]]
        .dropna(subset=["projid", "individualID"])
        .drop_duplicates(subset=["projid"])
    )
    out = longitudinal_df.copy()
    out = out.merge(mapping, on="projid", how="left")
    return out


def reshape_lipidomics_matrix(lipidomics_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert raw lipidomics matrix into sample-by-lipid format.

    Original table shape is lipid rows x sample columns with metadata columns
    ``label``, ``Mode``, and ``lipid class``. Output shape is sample rows x lipid
    columns with standardized names:

    ``{lipid_label}_{mode}_{lipid_class}``
    """
    required = {"label", "Mode", "lipid class"}
    missing = required.difference(lipidomics_df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"Missing required lipidomics columns: {missing_str}")

    metadata_cols = ["label", "Mode", "lipid class"]
    sample_cols = [c for c in lipidomics_df.columns if c not in metadata_cols]

    wide = lipidomics_df.set_index(metadata_cols)[sample_cols].T
    wide.index.name = "label"
    wide.columns = [
        f"{str(lipid)}_{str(mode)}_{str(lipid_class)}"
        for lipid, mode, lipid_class in wide.columns
    ]
    # Preserve naming correction from original notebook.
    wide.columns = [c.replace("postive", "positive") for c in wide.columns]

    out = wide.reset_index()
    lipid_cols = [c for c in out.columns if c != "label"]
    out[lipid_cols] = out[lipid_cols].apply(pd.to_numeric, errors="coerce")
    return out


def extract_individual_id_from_label(label: object) -> str | float:
    """
    Parse individual ID from the sample ``label`` field.

    The original notebooks expected labels like ``b1.127C_R8316516``.
    """
    if pd.isna(label):
        return np.nan
    text = str(label)
    if "_" not in text:
        return np.nan
    token = text.split("_")[-1]
    if token.startswith("R"):
        return token
    return np.nan


def add_individual_id_column(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Add ``individualID`` parsed from the sample label if not present."""
    out = sample_df.copy()
    if "individualID" not in out.columns:
        out["individualID"] = out["label"].apply(extract_individual_id_from_label)
    return out


def normalize_by_internal_standards(
    sample_df: pd.DataFrame,
    negative_ref_col: str = "PC(18:1D7_15:0)_negative_PC",
    positive_ref_col: str = "PC(18:1D7_15:0)_positive_PC",
) -> pd.DataFrame:
    """
    Normalize lipid intensities by mode-specific internal standards.

    This mirrors the mode split in the original notebook:
    - ``*_negative_*`` columns are divided by the negative internal standard.
    - ``*_positive_*`` columns are divided by the positive internal standard.
    """
    out = sample_df.copy()
    if negative_ref_col not in out.columns or positive_ref_col not in out.columns:
        return out

    negative_ref = out[negative_ref_col].replace(0, np.nan)
    positive_ref = out[positive_ref_col].replace(0, np.nan)

    for col in out.columns:
        if col in {"label", "individualID"}:
            continue
        if not is_lipid_column(col):
            continue

        parts = col.split("_")
        if len(parts) < 3:
            continue

        mode = parts[-2]
        if mode == "negative":
            out[col] = out[col] / negative_ref
        elif mode == "positive":
            out[col] = out[col] / positive_ref

    return out


def merge_lipid_id_map(sample_df: pd.DataFrame, lipid_individual_map_df: pd.DataFrame) -> pd.DataFrame:
    """Merge ``specimenID -> individualID`` mapping into the sample dataframe."""
    map_cols = ["specimenID", "individualID"]
    available_cols = [c for c in map_cols if c in lipid_individual_map_df.columns]
    if available_cols != map_cols:
        return sample_df.copy()

    out = sample_df.merge(
        lipid_individual_map_df[map_cols].drop_duplicates("specimenID"),
        left_on="label",
        right_on="specimenID",
        how="left",
        suffixes=("", "_mapped"),
    )
    out["individualID"] = out["individualID"].fillna(out["individualID_mapped"])
    out = out.drop(columns=[c for c in ["specimenID", "individualID_mapped"] if c in out.columns])
    return out


def attach_covariates(
    sample_df: pd.DataFrame,
    longitudinal_df: pd.DataFrame,
    covariate_columns: Sequence[str] = DEFAULT_COVARIATE_COLUMNS,
) -> pd.DataFrame:
    """Attach covariates from longitudinal metadata using ``individualID``."""
    needed = ["individualID"] + [c for c in covariate_columns if c in longitudinal_df.columns]
    long_sub = longitudinal_df[needed].copy()
    # Keep first non-null row per individual after sorting by original row order.
    long_sub = long_sub.drop_duplicates(subset=["individualID"], keep="first")

    out = sample_df.merge(long_sub, on="individualID", how="left")
    out = out.rename(columns={k: v for k, v in RENAME_COVARIATES.items() if k in out.columns})
    return out


def is_lipid_column(column_name: str) -> bool:
    """Return True when a column follows the lipid naming convention."""
    if "_" not in column_name:
        return False
    suffix = column_name.split("_")[-1]
    return suffix in LIPID_CATEGORY_SUFFIXES


def get_lipid_columns(df: pd.DataFrame) -> list[str]:
    """Get lipid feature columns from a dataframe."""
    return [col for col in df.columns if is_lipid_column(col)]


def drop_optional_columns(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    """Drop optional columns if they exist."""
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out = out.drop(columns=col)
    return out


def prepare_analysis_dataset(
    lipidomics_raw_df: pd.DataFrame,
    longitudinal_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    lipid_individual_map_df: pd.DataFrame,
    normalize_internal_standards: bool = True,
) -> pd.DataFrame:
    """
    Build the analysis-ready sample-by-lipid dataset from raw inputs.

    Returns a dataframe equivalent in role to the original
    ``Normalized_Formatted_Lipidomics.csv``/``Final_Formated_Lipidomics.csv``.
    """
    longitudinal_with_ids = add_individual_ids_to_longitudinal(longitudinal_df, clinical_df)
    sample_df = reshape_lipidomics_matrix(lipidomics_raw_df)
    sample_df = add_individual_id_column(sample_df)

    if normalize_internal_standards:
        sample_df = normalize_by_internal_standards(sample_df)

    sample_df = merge_lipid_id_map(sample_df, lipid_individual_map_df)
    sample_df = attach_covariates(sample_df, longitudinal_with_ids)
    sample_df = drop_optional_columns(sample_df, columns=["Net_SI", "Sex"])
    return sample_df

