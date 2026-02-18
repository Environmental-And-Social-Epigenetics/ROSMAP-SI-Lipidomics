"""Statistical helper functions shared across scripts and notebooks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import shapiro
from sklearn.neighbors import LocalOutlierFactor
from statsmodels.stats.multitest import fdrcorrection

from src.data_utils import get_lipid_columns


@dataclass(frozen=True)
class RegressionSpec:
    """Configuration for one per-lipid regression pass."""

    predictors: tuple[str, ...] = ("SI_avg", "niareagansc", "age_death")
    primary_predictor: str = "SI_avg"
    min_n: int = 20


def _q(column_name: str) -> str:
    """Safely quote dataframe column names for patsy/statsmodels formulas."""
    escaped = column_name.replace("\\", "\\\\").replace('"', '\\"')
    return f'Q("{escaped}")'


def build_formula(outcome: str, predictors: Sequence[str]) -> str:
    """Build a statsmodels-safe linear regression formula."""
    rhs = " + ".join(_q(p) for p in predictors)
    return f"{_q(outcome)} ~ {rhs}"


def clean_numeric_frame(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    """Keep requested columns and drop rows with missing or infinite values."""
    out = df[list(columns)].copy()
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(axis=0, how="any")
    return out


def run_lof_outlier_detection(
    df: pd.DataFrame,
    lipid_columns: Sequence[str] | None = None,
    n_neighbors: int = 100,
    strong_threshold: float = -1.5,
) -> pd.DataFrame:
    """
    Compute Local Outlier Factor scores for sample-level outlier checking.

    Returns one row per sample with:
    - ``lof_prediction`` (+1 inlier, -1 outlier)
    - ``lof_score`` (negative outlier factor; more negative is more outlier-like)
    - ``is_strong_outlier`` (score < ``strong_threshold``)
    """
    if lipid_columns is None:
        lipid_columns = get_lipid_columns(df)

    x = clean_numeric_frame(df, lipid_columns)
    if x.empty:
        return pd.DataFrame(columns=["index", "lof_prediction", "lof_score", "is_strong_outlier"])

    n_neighbors = max(2, min(n_neighbors, len(x) - 1))
    lof = LocalOutlierFactor(n_neighbors=n_neighbors, algorithm="brute")
    predictions = lof.fit_predict(x)
    scores = lof.negative_outlier_factor_

    result = pd.DataFrame(
        {
            "index": x.index,
            "lof_prediction": predictions,
            "lof_score": scores,
            "is_strong_outlier": scores < strong_threshold,
        }
    )
    return result


def run_shapiro_normality_tests(
    df: pd.DataFrame,
    lipid_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Run Shapiro-Wilk normality tests for each lipid column."""
    if lipid_columns is None:
        lipid_columns = get_lipid_columns(df)

    rows: list[dict[str, float | str | int]] = []
    for lipid in lipid_columns:
        values = pd.to_numeric(df[lipid], errors="coerce").dropna()
        if len(values) < 3:
            continue
        stat, p_value = shapiro(values)
        rows.append(
            {
                "lipid": lipid,
                "shapiro_stat": stat,
                "p_value": p_value,
                "n_samples": len(values),
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = add_fdr_column(out, p_value_column="p_value", output_column="fdr_p_value")
    return out


def add_fdr_column(
    df: pd.DataFrame,
    p_value_column: str = "p_value",
    output_column: str = "fdr_p_value",
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Add Benjamini-Hochberg FDR-adjusted p-values."""
    out = df.copy()
    if out.empty or p_value_column not in out.columns:
        out[output_column] = np.nan
        out[f"{output_column}_reject"] = False
        return out

    rejected, corrected = fdrcorrection(out[p_value_column].values, alpha=alpha, method="indep")
    out[output_column] = corrected
    out[f"{output_column}_reject"] = rejected
    return out


def run_per_lipid_regression(
    df: pd.DataFrame,
    lipid_columns: Sequence[str] | None = None,
    spec: RegressionSpec = RegressionSpec(),
) -> pd.DataFrame:
    """Fit OLS for each lipid and return coefficient/p-value summaries."""
    if lipid_columns is None:
        lipid_columns = get_lipid_columns(df)

    results: list[dict[str, float | str | int]] = []
    needed_predictors = list(spec.predictors)

    for lipid in lipid_columns:
        needed_cols = [lipid] + needed_predictors
        if any(col not in df.columns for col in needed_cols):
            continue

        model_df = clean_numeric_frame(df, needed_cols)
        if len(model_df) < spec.min_n:
            continue
        if model_df[lipid].nunique(dropna=True) <= 1:
            continue

        formula = build_formula(lipid, spec.predictors)
        try:
            model = smf.ols(formula=formula, data=model_df).fit()
        except Exception:
            continue

        primary_key = _q(spec.primary_predictor)
        if primary_key not in model.params.index:
            continue

        results.append(
            {
                "lipid": lipid,
                "n_samples": int(model.nobs),
                "r_squared": float(model.rsquared),
                "coef_primary": float(model.params[primary_key]),
                "p_value_primary": float(model.pvalues[primary_key]),
            }
        )

    out = pd.DataFrame(results)
    if not out.empty:
        out = add_fdr_column(out, p_value_column="p_value_primary", output_column="fdr_p_value")
        out = out.sort_values("p_value_primary", ascending=True).reset_index(drop=True)
    return out


def run_ancova_sex_interaction(
    df: pd.DataFrame,
    lipid_columns: Sequence[str] | None = None,
    main_predictor: str = "SI_avg",
    sex_column: str = "msex",
    covariates: Sequence[str] = ("niareagansc", "age_death"),
    min_n: int = 20,
) -> pd.DataFrame:
    """
    Run ANCOVA-style interaction models across lipids.

    Model form:
        lipid ~ SI_avg * msex + niareagansc + age_death

    The interaction term tests whether SI-lipid associations differ by sex.
    """
    if lipid_columns is None:
        lipid_columns = get_lipid_columns(df)

    required_shared = [main_predictor, sex_column] + list(covariates)
    if any(col not in df.columns for col in required_shared):
        return pd.DataFrame(
            columns=[
                "lipid",
                "n_samples",
                "r_squared",
                "coef_SI_avg",
                "p_SI_avg",
                "coef_msex",
                "p_msex",
                "coef_interaction",
                "p_interaction",
                "fdr_p_interaction",
                "fdr_p_interaction_reject",
            ]
        )

    main_key = _q(main_predictor)
    sex_key = _q(sex_column)
    interaction_key = f"{main_key}:{sex_key}"

    results: list[dict[str, float | str | int]] = []
    for lipid in lipid_columns:
        needed_cols = [lipid] + required_shared
        if any(col not in df.columns for col in needed_cols):
            continue

        model_df = clean_numeric_frame(df, needed_cols)
        if len(model_df) < min_n:
            continue
        if model_df[lipid].nunique(dropna=True) <= 1:
            continue

        covariate_terms = " + ".join(_q(c) for c in covariates)
        formula = f"{_q(lipid)} ~ {main_key} * {sex_key} + {covariate_terms}"
        try:
            model = smf.ols(formula=formula, data=model_df).fit()
        except Exception:
            continue

        if any(key not in model.params.index for key in (main_key, sex_key, interaction_key)):
            continue

        results.append(
            {
                "lipid": lipid,
                "n_samples": int(model.nobs),
                "r_squared": float(model.rsquared),
                "coef_SI_avg": float(model.params[main_key]),
                "p_SI_avg": float(model.pvalues[main_key]),
                "coef_msex": float(model.params[sex_key]),
                "p_msex": float(model.pvalues[sex_key]),
                "coef_interaction": float(model.params[interaction_key]),
                "p_interaction": float(model.pvalues[interaction_key]),
            }
        )

    out = pd.DataFrame(results)
    if not out.empty:
        out = add_fdr_column(
            out,
            p_value_column="p_interaction",
            output_column="fdr_p_interaction",
        )
        out = out.sort_values("p_interaction", ascending=True).reset_index(drop=True)
    return out


def split_by_sex(df: pd.DataFrame, sex_column: str = "msex") -> dict[str, pd.DataFrame]:
    """Return male/female/full cohorts keyed by descriptive labels."""
    cohorts = {"all": df.copy()}
    if sex_column in df.columns:
        cohorts["male"] = df[df[sex_column] == 1].copy()
        cohorts["female"] = df[df[sex_column] == 0].copy()
    return cohorts


def filter_no_ad(df: pd.DataFrame, reagan_column: str = "niareagansc", threshold: float = 2) -> pd.DataFrame:
    """Filter to samples with no/low AD pathology per original notebook logic."""
    if reagan_column not in df.columns:
        return df.copy()
    return df[pd.to_numeric(df[reagan_column], errors="coerce") > threshold].copy()


def infer_lipid_categories(lipid_columns: Sequence[str]) -> list[str]:
    """Infer lipid categories from suffixes in lipid column names."""
    categories = sorted({c.split("_")[-1] for c in lipid_columns if "_" in c})
    return categories


def zscore_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Column-wise z-score transform for selected numeric columns."""
    out = df.copy()
    for col in columns:
        series = pd.to_numeric(out[col], errors="coerce")
        sd = series.std(ddof=0)
        if pd.isna(sd) or sd == 0:
            out[col] = np.nan
            continue
        out[col] = (series - series.mean()) / sd
    return out


def compute_category_means(
    df: pd.DataFrame,
    lipid_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Add one averaged column per lipid category.

    Example output columns: ``catmean_Cer``, ``catmean_TG``, etc.
    """
    out = df.copy()
    if lipid_columns is None:
        lipid_columns = get_lipid_columns(df)

    category_map: dict[str, list[str]] = {}
    for col in lipid_columns:
        cat = col.split("_")[-1]
        category_map.setdefault(cat, []).append(col)

    for category, cols in sorted(category_map.items()):
        out[f"catmean_{category}"] = out[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)

    return out

