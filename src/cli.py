"""
Command-line interface and main entry point for AREA.

Usage
-----
    # All parameters on the command line:
    python run_area.py -bf bools.csv -rf ranks.csv -jc sample_id -od results/

    # All parameters from a config file:
    python run_area.py --config area_config.yaml

    # Config file with command-line overrides (CLI wins):
    python run_area.py --config area_config.yaml --verbose --threads 8
"""

import argparse
import os
import time

import pandas as pd

from .config import load_config
from .correction import compute_adjusted_pvalues
from .planning import build_run_plan
from .runner import compute_pvalues


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    """Build and parse the command-line interface.

    If ``--config`` is supplied, its values are loaded first.  Any
    explicit CLI flags override the corresponding config values.
    """
    parser = argparse.ArgumentParser(
        description="Run AREA (Association by Rank Enrichment Analysis)"
    )

    parser.add_argument("-c", "--config", default=None,
                        help="Path to a YAML config file")

    parser.add_argument("-bf", "--boolean-file", default=None,
                        help="CSV of binary attributes (samples x attributes)")
    parser.add_argument("-rf", "--rank-file", default=None,
                        help="CSV of continuous values (samples x features)")
    parser.add_argument("-jc", "--join-column", default=None,
                        help="Column shared between both input files")
    parser.add_argument("-od", "--out-dir", default=None,
                        help="Directory for all output files")
    parser.add_argument("-t", "--threads", type=int, default=None,
                        help="Number of parallel threads (default: 4)")

    parser.add_argument("--keep-rank-columns", default=None,
                        help="Text file listing rank columns to include")
    parser.add_argument("--keep-bool-columns", default=None,
                        help="Text file listing bool columns to include")
    parser.add_argument("--keep-samples", default=None,
                        help="Text file listing sample IDs to include")

    parser.add_argument("--gpu", action="store_true", default=None,
                        help="Use GPU acceleration via CuPy")
    parser.add_argument("--verbose", action="store_true", default=None,
                        help="Print diagnostic output")

    args = parser.parse_args(argv)

    # -- Merge config file with CLI flags -------------------------------------
    # Strategy: start with defaults, layer on config, layer on CLI.
    defaults = {
        "boolean_file": None,
        "rank_file": None,
        "join_column": None,
        "out_dir": None,
        "threads": 4,
        "keep_rank_columns": None,
        "keep_bool_columns": None,
        "keep_samples": None,
        "gpu": False,
        "verbose": False,
    }

    # Layer 1: config file
    if args.config is not None:
        config_values = load_config(args.config)
        defaults.update(config_values)

    # Layer 2: explicit CLI flags override config.
    # argparse sets unspecified flags to None (or None for store_true
    # with default=None), so we only override when the user actually
    # passed a flag.
    cli_dict = vars(args)
    for key in defaults:
        cli_value = cli_dict.get(key)
        if cli_value is not None:
            defaults[key] = cli_value

    # Validate required fields
    required = ["boolean_file", "rank_file", "join_column", "out_dir"]
    missing = [k for k in required if not defaults[k]]
    if missing:
        parser.error(
            f"Missing required parameter(s): {', '.join(missing)}. "
            "Provide them via --config or the command line."
        )

    # Return as a Namespace so the rest of the code can use args.X
    return argparse.Namespace(**defaults)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    plan_file = f"area_scores_{timestamp}"

    # -- Array backend --------------------------------------------------------
    backend_name = "cupy (GPU)" if args.gpu else "numpy (CPU)"
    print(f"Array backend: {backend_name}")

    # -- Read inputs ----------------------------------------------------------
    bool_df = pd.read_csv(args.boolean_file, index_col=0).reset_index()
    rank_df = pd.read_csv(args.rank_file, index_col=0).reset_index()

    # -- Sample keep-list -----------------------------------------------------
    if args.keep_samples is not None:
        keep_samples = pd.read_csv(args.keep_samples, header=None)[0].tolist()
    else:
        keep_samples = []

    # -- Step 1: build plan ---------------------------------------------------
    print("=== Building run plan ===")
    build_run_plan(
        out_dir=args.out_dir,
        plan_file=plan_file,
        join_column=args.join_column,
        rank_df=rank_df,
        bool_df=bool_df,
        keep_bools_file=args.keep_bool_columns,
        keep_ranks_file=args.keep_rank_columns,
        verbose=args.verbose,
    )

    # -- Step 2: compute p-values ---------------------------------------------
    print("=== Computing enrichment scores ===")
    compute_pvalues(
        out_dir=args.out_dir,
        plan_file=plan_file,
        join_column=args.join_column,
        rank_df=rank_df,
        bool_df=bool_df,
        keep_samples=keep_samples,
        use_gpu=args.gpu,
        n_threads=args.threads,
        verbose=args.verbose,
    )

    # -- Step 3: adjust p-values ----------------------------------------------
    print("=== Adjusting p-values ===")
    compute_adjusted_pvalues(out_dir=args.out_dir, plan_file=plan_file)

    # -- Cleanup intermediate files -------------------------------------------
    if not args.verbose:
        plan_path = os.path.join(args.out_dir, plan_file)
        raw_path = os.path.join(args.out_dir, plan_file + ".raw_pvalues.csv")
        for path in (plan_path, raw_path):
            if os.path.isfile(path):
                os.remove(path)
    else:
        print("Verbose mode: keeping intermediate files "
              f"({plan_file}, {plan_file}.raw_pvalues.csv)")

    print(f"AREA analysis complete. Output prefix: {plan_file}")


if __name__ == "__main__":
    main()
