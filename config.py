"""Central configuration for repository paths and common file names."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

DATA_DIR = PROJECT_ROOT / "data"
DATA_RAW_DIR = DATA_DIR / "raw"
DATA_PROCESSED_DIR = DATA_DIR / "processed"

RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
SRC_DIR = PROJECT_ROOT / "src"

# Canonical input file names copied from the original analysis workspace.
RAW_LIPIDOMICS_FILENAME = "ROSMAP_Brain_SERRF_normalization_internal.csv"
META_LONGITUDINAL_FILENAME = "dataset_652_basic_03-23-2022.csv"
META_CLINICAL_FILENAME = "ROSMAP_clinical.csv"
LIPID_ID_MAP_FILENAME = "Lip_Ind_Map.csv"
ORIGINAL_FINAL_FORMATTED_FILENAME = "Final_Formated_Lipidomics.csv"
ORIGINAL_NORMALIZED_FORMATTED_FILENAME = "Normalized_Formatted_Lipidomics.csv"

# Canonical processed output names used by this rebuilt repository.
FINAL_FORMATTED_FILENAME = "Final_Formatted_Lipidomics.csv"
NORMALIZED_FORMATTED_FILENAME = "Normalized_Formatted_Lipidomics.csv"


def ensure_project_dirs() -> None:
    """Create the expected local folder tree if it does not exist."""
    for directory in (
        DATA_RAW_DIR,
        DATA_PROCESSED_DIR,
        FIGURES_DIR,
        TABLES_DIR,
        NOTEBOOKS_DIR,
        SCRIPTS_DIR,
        SRC_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)
