"""
Run-plan construction: determine which (rank_column, bool_column) pairs
to analyse, applying keep-file filters and dropping constant columns.
"""

import os

import pandas as pd


def build_run_plan(out_dir, plan_file, join_column, rank_df, bool_df,
                   keep_bools_file=None, keep_ranks_file=None,
                   verbose=False):
    """Determine which (rank_column, bool_column) pairs to analyse.

    Writes a CSV plan file to ``out_dir/plan_file`` listing every valid
    pair with ``plan = run_area``.  Pairs involving constant columns or
    columns excluded by keep-files are omitted entirely.

    Parameters
    ----------
    out_dir : str
        Output directory (must exist).
    plan_file : str
        Filename for the plan CSV.
    join_column : str
        Column shared between rank and bool DataFrames.
    rank_df : DataFrame
        Continuous-value data (samples x features).
    bool_df : DataFrame
        Binary-attribute data (samples x attributes).
    keep_bools_file : str or None
        Path to a single-column text file listing bool columns to retain.
    keep_ranks_file : str or None
        Path to a single-column text file listing rank columns to retain.
    verbose : bool
        Print diagnostic counts.
    """
    rank_columns = [c for c in rank_df.columns if c != join_column]
    bool_columns = [c for c in bool_df.columns if c != join_column]

    if verbose:
        print(f"Initial rank columns: {len(rank_columns)}")
        print(f"Initial bool columns: {len(bool_columns)}")

    # Apply keep-file filters
    if keep_bools_file is not None:
        allowed = set(pd.read_csv(keep_bools_file, header=None)[0])
        bool_columns = [b for b in bool_columns if b in allowed]
        if verbose:
            print(f"After keep_bools filter: {len(bool_columns)} bool columns")

    if keep_ranks_file is not None:
        allowed = set(pd.read_csv(keep_ranks_file, header=None)[0])
        rank_columns = [r for r in rank_columns if r in allowed]
        if verbose:
            print(f"After keep_ranks filter: {len(rank_columns)} rank columns")

    # Drop constant columns
    rank_columns = [c for c in rank_columns
                    if rank_df[c].nunique(dropna=False) > 1]
    bool_columns = [c for c in bool_columns
                    if bool_df[c].nunique(dropna=False) > 1]

    if verbose:
        print(f"After removing constant columns: "
              f"{len(rank_columns)} rank, {len(bool_columns)} bool")
        print(f"Total pair count: {len(rank_columns) * len(bool_columns):,}")

    # Stream plan to disk
    plan_path = os.path.join(out_dir, plan_file)
    with open(plan_path, "w") as fh:
        fh.write("rank_column,bool_column,plan\n")
        for rank_col in rank_columns:
            chunk = pd.DataFrame({
                "rank_column": [rank_col] * len(bool_columns),
                "bool_column": bool_columns,
                "plan": "run_area",
            })
            chunk.to_csv(fh, header=False, index=False)
