"""
Orchestration: merge data, run enrichment per boolean column, and
coordinate parallel execution via threads.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from .backend import load_array_backend
from .enrichment import (
    compute_enrichment_score,
    compute_nes_pvalue,
    permute_enrichment_scores,
)


# ---------------------------------------------------------------------------
# Single boolean-column worker
# ---------------------------------------------------------------------------

def _run_single_bool_column(bool_col, rank_columns, join_column,
                            rank_df, bool_df, keep_samples,
                            use_gpu, n_permutations=1000,
                            verbose=False):
    """Run AREA for one boolean attribute against all its rank columns.

    Parameters
    ----------
    bool_col : str
        Name of the boolean attribute column.
    rank_columns : list[str]
        Rank columns to test against this attribute.
    join_column : str
        Column shared between both DataFrames.
    rank_df : DataFrame
        Continuous-value data.
    bool_df : DataFrame
        Binary-attribute data.
    keep_samples : list
        Sample IDs to retain (empty list means keep all).
    use_gpu : bool
        Use CuPy instead of NumPy.
    n_permutations : int
        Number of permutations for the null distribution.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    DataFrame
        Columns: [bool_column, rank_column, nes, pvalue]
    """
    xp = load_array_backend(use_gpu)

    if verbose:
        print(f"Processing bool column '{bool_col}' "
              f"against {len(rank_columns)} rank columns")

    # Subset and merge
    cols_bool = [bool_col, join_column]
    cols_rank = rank_columns + [join_column]
    merged = bool_df[cols_bool].merge(rank_df[cols_rank],
                                      on=join_column, how="inner")

    if len(keep_samples) > 0:
        merged = merged[merged[join_column].isin(keep_samples)]

    # Shuffle row order so tied-rank zeros are in random order
    merged = merged.sample(frac=1).reset_index(drop=True)

    # Build null distribution once per boolean attribute
    binary_vector = merged[bool_col].tolist()
    null_scores = permute_enrichment_scores(
        binary_vector, n_permutations=n_permutations, xp=xp, verbose=verbose,
    )

    # Score each rank column
    results = []
    for rank_col in rank_columns:
        sorted_df = merged.sort_values(rank_col)
        binary_sorted = sorted_df[bool_col].tolist()

        observed_es, *_ = compute_enrichment_score(binary_sorted, xp=xp)
        nes, pvalue = compute_nes_pvalue(observed_es, null_scores,
                                         use_gpu=use_gpu, xp=xp)
        results.append({
            "bool_column": bool_col,
            "rank_column": rank_col,
            "nes": nes,
            "pvalue": pvalue,
        })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Thread-pool dispatcher
# ---------------------------------------------------------------------------

def compute_pvalues(out_dir, plan_file, join_column, rank_df, bool_df,
                    keep_samples, use_gpu=False, n_threads=4,
                    verbose=False):
    """Compute NES + p-values for every planned pair.

    Uses threads instead of processes to avoid pickling.  NumPy and CuPy
    release the GIL during C-level computation, so threads provide real
    parallelism for the array-heavy workload.

    Results are written to ``<out_dir>/<plan_file>.raw_pvalues.csv``.
    """
    plan_path = os.path.join(out_dir, plan_file)
    plan_df = pd.read_csv(plan_path)
    plan_df = plan_df[plan_df["plan"] == "run_area"]

    bool_columns = plan_df["bool_column"].unique()
    print(f"Computing p-values for {len(bool_columns)} boolean attributes …")

    results = []
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {
            executor.submit(
                _run_single_bool_column,
                bool_col=bool_col,
                rank_columns=plan_df.loc[
                    plan_df["bool_column"] == bool_col, "rank_column"
                ].tolist(),
                join_column=join_column,
                rank_df=rank_df,
                bool_df=bool_df,
                keep_samples=keep_samples,
                use_gpu=use_gpu,
                verbose=verbose,
            ): bool_col
            for bool_col in bool_columns
        }
        for future in as_completed(futures):
            results.append(future.result())

    combined = pd.concat(results, ignore_index=True)
    raw_path = os.path.join(out_dir, plan_file + ".raw_pvalues.csv")
    combined.to_csv(raw_path, index=False)
    print("Done computing p-values.")
