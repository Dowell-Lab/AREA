"""
Multiple-testing p-value correction (Bonferroni, Holm, BH, BY).
"""

import os

import pandas as pd
from statsmodels.stats.multitest import multipletests


def _adjust(pvalues, method):
    """Apply a single correction method via ``multipletests``."""
    return multipletests(pvalues, method=method)[1]


def compute_adjusted_pvalues(out_dir, plan_file):
    """Read raw p-values, add multiple-testing corrections, write final CSV.

    Reads  ``<out_dir>/<plan_file>.raw_pvalues.csv``
    Writes ``<out_dir>/<plan_file>.adjusted_pvalues.csv``

    Rows whose p-value could not be parsed as numeric are dropped.
    """
    raw_path = os.path.join(out_dir, plan_file + ".raw_pvalues.csv")
    df = pd.read_csv(raw_path)

    # Coerce to numeric; rows that fail become NaN and are excluded
    df["pvalue_numeric"] = pd.to_numeric(df["pvalue"], errors="coerce")
    numeric_df = df.loc[df["pvalue_numeric"].notna()].copy()
    numeric_df["pvalue"] = numeric_df["pvalue_numeric"]

    pvals = numeric_df["pvalue"]
    numeric_df["pvalue_bonferroni"] = _adjust(pvals, "bonferroni")
    numeric_df["pvalue_holm"] = _adjust(pvals, "holm")
    numeric_df["pvalue_bh"] = _adjust(pvals, "fdr_bh")
    numeric_df["pvalue_by"] = _adjust(pvals, "fdr_by")

    numeric_df = numeric_df.drop(columns=["pvalue_numeric"])
    numeric_df = numeric_df.sort_values("pvalue_bh")

    out_path = os.path.join(out_dir, plan_file + ".adjusted_pvalues.csv")
    numeric_df.to_csv(out_path, index=False)
    print(f"Adjusted p-values written to {out_path}")
