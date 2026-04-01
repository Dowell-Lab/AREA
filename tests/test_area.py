"""
Unit and integration tests for the AREA package.

Tests use real data from ``AREA/testdata/``:
    - value_expression.csv — rank file   (samples × genes)
    - comorbid_file.csv    — boolean file (samples × comorbidities)

Run from the project root:
    python run_tests.py
    python run_tests.py -v    # verbose output
"""

import os
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_HERE)          # AREA/
_TESTDATA = os.path.join(_PROJECT_ROOT, "testdata")

RANK_FILE = os.path.join(_TESTDATA, "value_expression.csv")
BOOL_FILE = os.path.join(_TESTDATA, "comorbid_file.csv")


def _data_available():
    return os.path.isfile(RANK_FILE) and os.path.isfile(BOOL_FILE)


# The join column shared between both test CSV files.
# This must be specified explicitly — AREA never guesses the join column.
JOIN_COLUMN = "Patient"


# ===================================================================
# 1) Backend tests
# ===================================================================

class TestBackend(unittest.TestCase):
    """Verify array-backend helpers."""

    def test_load_cpu(self):
        from area.backend import load_array_backend
        xp = load_array_backend(use_gpu=False)
        self.assertIs(xp, np)

    def test_to_numpy_cpu(self):
        from area.backend import to_numpy
        val = np.float64(3.14)
        self.assertEqual(to_numpy(val, use_gpu=False), val)

    def test_trapz_cpu(self):
        from area.backend import trapz
        y = np.array([0.0, 0.5, 1.0])
        result = trapz(y, np)
        # Trapezoidal rule: (0+0.5)/2 * 1 + (0.5+1.0)/2 * 1 = 1.0
        self.assertAlmostEqual(float(result), 1.0, places=5)


# ===================================================================
# 2) Config tests
# ===================================================================

class TestConfig(unittest.TestCase):
    """YAML config loading."""

    def test_load_valid_config(self):
        from area.config import load_config
        with tempfile.NamedTemporaryFile("w", suffix=".yaml",
                                         delete=False) as f:
            f.write(
                "boolean_file: bools.csv\n"
                "rank_file: ranks.csv\n"
                "join_column: Patient\n"
                "out_dir: /tmp/out\n"
                "threads: 8\n"
                "verbose: true\n"
            )
            f.flush()
            cfg = load_config(f.name)
        os.unlink(f.name)

        self.assertEqual(cfg["boolean_file"], "bools.csv")
        self.assertEqual(cfg["rank_file"], "ranks.csv")
        self.assertEqual(cfg["join_column"], "Patient")
        self.assertEqual(cfg["threads"], 8)
        self.assertTrue(cfg["verbose"])

    def test_unknown_key_raises(self):
        from area.config import load_config
        with tempfile.NamedTemporaryFile("w", suffix=".yaml",
                                         delete=False) as f:
            f.write("unknown_key: 123\n")
            f.flush()
            with self.assertRaises(ValueError):
                load_config(f.name)
        os.unlink(f.name)


# ===================================================================
# 3) Enrichment math tests
# ===================================================================

class TestEnrichment(unittest.TestCase):
    """Core enrichment-score and NES maths."""

    def test_enrichment_score_perfect(self):
        """All hits at the front → positive enrichment score."""
        from area.enrichment import compute_enrichment_score
        binary = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        es, norm, trend, cumul = compute_enrichment_score(binary, xp=np)
        self.assertGreater(float(es), 0,
                           "Hits at front should yield positive ES")

    def test_enrichment_score_perfect_tail(self):
        """All hits at the tail → negative enrichment score."""
        from area.enrichment import compute_enrichment_score
        binary = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
        es, *_ = compute_enrichment_score(binary, xp=np)
        self.assertLess(float(es), 0,
                        "Hits at tail should yield negative ES")

    def test_enrichment_score_uniform(self):
        """Evenly spaced hits → ES near zero (within ±0.2)."""
        from area.enrichment import compute_enrichment_score
        binary = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        es, *_ = compute_enrichment_score(binary, xp=np)
        self.assertAlmostEqual(float(es), 0.0, delta=0.2,
                               msg="Uniform hits should give ~0 ES")

    def test_permutation_length(self):
        """Permutation null returns the right number of scores."""
        from area.enrichment import permute_enrichment_scores
        binary = [1, 0, 1, 0, 0]
        null = permute_enrichment_scores(binary, n_permutations=50,
                                         seed=42, xp=np)
        self.assertEqual(len(null), 50)

    def test_nes_pvalue_range(self):
        """NES p-value should be in [0, 1]."""
        from area.enrichment import (
            compute_enrichment_score,
            compute_nes_pvalue,
            permute_enrichment_scores,
        )
        binary = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        es, *_ = compute_enrichment_score(binary, xp=np)
        null = permute_enrichment_scores(binary, n_permutations=200,
                                         seed=42, xp=np)
        nes, pval = compute_nes_pvalue(es, null, use_gpu=False, xp=np)
        self.assertTrue(0 <= pval <= 1, f"p-value {pval} out of range")
        self.assertIsInstance(nes, float)


# ===================================================================
# 4) Planning tests (uses real data)
# ===================================================================

@unittest.skipUnless(_data_available(), "testdata/ not found — skipping")
class TestPlanning(unittest.TestCase):
    """Run-plan generation against the real test data."""

    @classmethod
    def setUpClass(cls):
        cls.rank_df = pd.read_csv(RANK_FILE, index_col=0).reset_index()
        cls.bool_df = pd.read_csv(BOOL_FILE, index_col=0).reset_index()
        cls.join_column = JOIN_COLUMN
        cls.tmpdir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_plan_file_created(self):
        from area.planning import build_run_plan
        plan_name = "test_plan"
        build_run_plan(
            out_dir=self.tmpdir,
            plan_file=plan_name,
            join_column=self.join_column,
            rank_df=self.rank_df,
            bool_df=self.bool_df,
            verbose=False,
        )
        plan_path = os.path.join(self.tmpdir, plan_name)
        self.assertTrue(os.path.isfile(plan_path),
                        f"Plan file not created at {plan_path}")

    def test_plan_has_rows(self):
        from area.planning import build_run_plan
        plan_name = "test_plan_rows"
        build_run_plan(
            out_dir=self.tmpdir,
            plan_file=plan_name,
            join_column=self.join_column,
            rank_df=self.rank_df,
            bool_df=self.bool_df,
        )
        plan_df = pd.read_csv(os.path.join(self.tmpdir, plan_name))
        self.assertGreater(len(plan_df), 0, "Plan should have rows")
        self.assertListEqual(
            sorted(plan_df.columns.tolist()),
            sorted(["rank_column", "bool_column", "plan"]),
        )
        self.assertTrue((plan_df["plan"] == "run_area").all())

    def test_plan_skips_constant_columns(self):
        """A rank column with all-same values should be excluded."""
        from area.planning import build_run_plan
        rank_modified = self.rank_df.copy()
        # Make one column constant
        const_col = [c for c in rank_modified.columns
                     if c != self.join_column][0]
        rank_modified[const_col] = 99.0

        plan_name = "test_plan_const"
        build_run_plan(
            out_dir=self.tmpdir,
            plan_file=plan_name,
            join_column=self.join_column,
            rank_df=rank_modified,
            bool_df=self.bool_df,
        )
        plan_df = pd.read_csv(os.path.join(self.tmpdir, plan_name))
        self.assertNotIn(const_col, plan_df["rank_column"].values,
                         "Constant column should be excluded from plan")


# ===================================================================
# 5) Runner tests (uses real data)
# ===================================================================

@unittest.skipUnless(_data_available(), "testdata/ not found — skipping")
class TestRunner(unittest.TestCase):
    """Thread-pool p-value computation with the real test files."""

    @classmethod
    def setUpClass(cls):
        cls.rank_df = pd.read_csv(RANK_FILE, index_col=0).reset_index()
        cls.bool_df = pd.read_csv(BOOL_FILE, index_col=0).reset_index()
        cls.join_column = JOIN_COLUMN
        cls.tmpdir = tempfile.mkdtemp()

        # Build a small plan first (limit to 2 bool cols × 2 rank cols)
        rank_cols = [c for c in cls.rank_df.columns
                     if c != cls.join_column][:2]
        bool_cols = [c for c in cls.bool_df.columns
                     if c != cls.join_column
                     and cls.bool_df[c].nunique(dropna=False) > 1][:2]

        cls.plan_file = "runner_test_plan"
        plan_path = os.path.join(cls.tmpdir, cls.plan_file)
        rows = []
        for rc in rank_cols:
            for bc in bool_cols:
                rows.append({"rank_column": rc, "bool_column": bc,
                             "plan": "run_area"})
        pd.DataFrame(rows).to_csv(plan_path, index=False)
        cls.rank_cols = rank_cols
        cls.bool_cols = bool_cols

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_compute_pvalues_creates_output(self):
        from area.runner import compute_pvalues
        compute_pvalues(
            out_dir=self.tmpdir,
            plan_file=self.plan_file,
            join_column=self.join_column,
            rank_df=self.rank_df,
            bool_df=self.bool_df,
            keep_samples=[],
            use_gpu=False,
            n_threads=2,
            verbose=False,
        )
        raw_path = os.path.join(
            self.tmpdir, self.plan_file + ".raw_pvalues.csv"
        )
        self.assertTrue(os.path.isfile(raw_path), "raw_pvalues CSV missing")

        results = pd.read_csv(raw_path)
        expected_rows = len(self.rank_cols) * len(self.bool_cols)
        self.assertEqual(len(results), expected_rows,
                         "Row count should equal rank × bool pairs")
        self.assertTrue(
            set(["bool_column", "rank_column", "nes", "pvalue"])
            .issubset(results.columns),
        )

    def test_pvalues_in_range(self):
        raw_path = os.path.join(
            self.tmpdir, self.plan_file + ".raw_pvalues.csv"
        )
        if not os.path.isfile(raw_path):
            self.skipTest("raw_pvalues not yet generated")
        results = pd.read_csv(raw_path)
        pvals = pd.to_numeric(results["pvalue"], errors="coerce").dropna()
        self.assertTrue((pvals >= 0).all() and (pvals <= 1).all(),
                        "All p-values should be in [0, 1]")


# ===================================================================
# 6) Correction tests (uses real data)
# ===================================================================

@unittest.skipUnless(_data_available(), "testdata/ not found — skipping")
class TestCorrection(unittest.TestCase):
    """Multiple-testing correction on real runner output."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.plan_file = "correction_test_plan"

        # Create a synthetic raw_pvalues file with realistic numbers
        raw_data = pd.DataFrame({
            "bool_column": ["A", "A", "B", "B"],
            "rank_column": ["g1", "g2", "g1", "g2"],
            "nes": [-1.2, 0.8, -0.5, 1.1],
            "pvalue": [0.001, 0.05, 0.30, 0.02],
        })
        raw_path = os.path.join(
            cls.tmpdir, cls.plan_file + ".raw_pvalues.csv"
        )
        raw_data.to_csv(raw_path, index=False)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_adjusted_file_created(self):
        from area.correction import compute_adjusted_pvalues
        compute_adjusted_pvalues(
            out_dir=self.tmpdir, plan_file=self.plan_file
        )
        adj_path = os.path.join(
            self.tmpdir, self.plan_file + ".adjusted_pvalues.csv"
        )
        self.assertTrue(os.path.isfile(adj_path))

    def test_adjusted_columns(self):
        adj_path = os.path.join(
            self.tmpdir, self.plan_file + ".adjusted_pvalues.csv"
        )
        if not os.path.isfile(adj_path):
            self.skipTest("adjusted_pvalues not yet generated")
        df = pd.read_csv(adj_path)
        for col in ["pvalue_bonferroni", "pvalue_holm",
                     "pvalue_bh", "pvalue_by"]:
            self.assertIn(col, df.columns, f"Missing column: {col}")

    def test_adjusted_pvalues_geq_raw(self):
        """Adjusted p-values should be ≥ raw p-values (corrections inflate)."""
        adj_path = os.path.join(
            self.tmpdir, self.plan_file + ".adjusted_pvalues.csv"
        )
        if not os.path.isfile(adj_path):
            self.skipTest("adjusted_pvalues not yet generated")
        df = pd.read_csv(adj_path)
        for col in ["pvalue_bonferroni", "pvalue_holm"]:
            self.assertTrue(
                (df[col] >= df["pvalue"] - 1e-12).all(),
                f"{col} should be ≥ raw p-value",
            )

    def test_sorted_by_bh(self):
        adj_path = os.path.join(
            self.tmpdir, self.plan_file + ".adjusted_pvalues.csv"
        )
        if not os.path.isfile(adj_path):
            self.skipTest("adjusted_pvalues not yet generated")
        df = pd.read_csv(adj_path)
        bh = df["pvalue_bh"].tolist()
        self.assertEqual(bh, sorted(bh),
                         "Output should be sorted by BH-adjusted p-value")


# ===================================================================
# 7) CLI / parse_args tests
# ===================================================================

class TestCLI(unittest.TestCase):
    """Argument parsing and config merging."""

    def test_parse_args_minimal(self):
        from area.cli import parse_args
        args = parse_args([
            "-bf", "bools.csv",
            "-rf", "ranks.csv",
            "-jc", "Patient",
            "-od", "/tmp/out",
        ])
        self.assertEqual(args.boolean_file, "bools.csv")
        self.assertEqual(args.rank_file, "ranks.csv")
        self.assertEqual(args.join_column, "Patient")
        self.assertEqual(args.out_dir, "/tmp/out")
        self.assertEqual(args.threads, 4)       # default
        self.assertFalse(args.gpu)               # default
        self.assertFalse(args.verbose)           # default

    def test_parse_args_missing_required(self):
        from area.cli import parse_args
        with self.assertRaises(SystemExit):
            parse_args(["-bf", "bools.csv"])     # missing rf, jc, od

    def test_parse_args_with_config(self):
        from area.cli import parse_args
        with tempfile.NamedTemporaryFile("w", suffix=".yaml",
                                         delete=False) as f:
            f.write(
                "boolean_file: from_config.csv\n"
                "rank_file: ranks.csv\n"
                "join_column: Patient\n"
                "out_dir: /tmp/cfg\n"
                "threads: 16\n"
            )
            f.flush()
            # CLI overrides config for boolean_file
            args = parse_args([
                "--config", f.name,
                "-bf", "from_cli.csv",
            ])
        os.unlink(f.name)
        self.assertEqual(args.boolean_file, "from_cli.csv",
                         "CLI should override config")
        self.assertEqual(args.threads, 16,
                         "Config value used when CLI absent")


# ===================================================================
# 8) End-to-end integration test (uses real data)
# ===================================================================

@unittest.skipUnless(_data_available(), "testdata/ not found — skipping")
class TestEndToEnd(unittest.TestCase):
    """Full pipeline via cli.main() using real test files."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_full_pipeline_verbose(self):
        """Run the entire AREA pipeline with --verbose (keeps intermediate files).

        Checks that:
        - The plan file is created with at least one pair
        - raw_pvalues CSV has correct columns and valid p-values
        - adjusted_pvalues CSV has correction columns and BH ordering
        - Intermediate files are preserved when --verbose is set
        """
        from area.cli import main
        verbose_dir = os.path.join(self.tmpdir, "verbose")
        os.makedirs(verbose_dir, exist_ok=True)

        main([
            "-bf", BOOL_FILE,
            "-rf", RANK_FILE,
            "-jc", self._detect_join_column(),
            "-od", verbose_dir,
            "-t", "2",
            "--verbose",
        ])

        # --- Plan file ---------------------------------------------------
        plan_files = [f for f in os.listdir(verbose_dir)
                      if not f.endswith(".csv")]
        self.assertTrue(len(plan_files) >= 1,
                        "At least one plan file should exist")
        plan_name = plan_files[0]
        plan_df = pd.read_csv(os.path.join(verbose_dir, plan_name))
        self.assertGreater(len(plan_df), 0, "Plan should have pairs")

        # --- Raw p-values ------------------------------------------------
        raw_path = os.path.join(verbose_dir, plan_name + ".raw_pvalues.csv")
        self.assertTrue(os.path.isfile(raw_path),
                        "raw_pvalues should be kept in verbose mode")
        raw_df = pd.read_csv(raw_path)
        self.assertGreater(len(raw_df), 0, "raw_pvalues should have rows")
        self.assertTrue(
            {"bool_column", "rank_column", "nes", "pvalue"}
            .issubset(raw_df.columns)
        )
        pvals = pd.to_numeric(raw_df["pvalue"], errors="coerce").dropna()
        self.assertTrue((pvals >= 0).all() and (pvals <= 1).all(),
                        "All raw p-values should be in [0, 1]")

        # --- Adjusted p-values -------------------------------------------
        adj_path = os.path.join(
            verbose_dir, plan_name + ".adjusted_pvalues.csv"
        )
        self.assertTrue(os.path.isfile(adj_path), "adjusted_pvalues missing")
        adj_df = pd.read_csv(adj_path)
        for col in ["pvalue_bonferroni", "pvalue_holm",
                     "pvalue_bh", "pvalue_by"]:
            self.assertIn(col, adj_df.columns, f"Missing column: {col}")
        # BH column should be sorted ascending
        bh = adj_df["pvalue_bh"].tolist()
        self.assertEqual(bh, sorted(bh), "Output should be sorted by BH")

    def test_cleanup_without_verbose(self):
        """Without --verbose, intermediate files should be deleted.

        Only the final adjusted_pvalues CSV should remain.
        """
        from area.cli import main
        clean_dir = os.path.join(self.tmpdir, "clean")
        os.makedirs(clean_dir, exist_ok=True)

        main([
            "-bf", BOOL_FILE,
            "-rf", RANK_FILE,
            "-jc", self._detect_join_column(),
            "-od", clean_dir,
            "-t", "2",
        ])

        remaining = os.listdir(clean_dir)
        # Only the adjusted_pvalues CSV should remain
        adj_files = [f for f in remaining
                     if f.endswith(".adjusted_pvalues.csv")]
        self.assertEqual(len(adj_files), 1,
                         "Exactly one adjusted_pvalues file expected")

        # Plan and raw files should be gone
        plan_files = [f for f in remaining if not f.endswith(".csv")]
        raw_files = [f for f in remaining
                     if f.endswith(".raw_pvalues.csv")]
        self.assertEqual(len(plan_files), 0,
                         "Plan file should be deleted without --verbose")
        self.assertEqual(len(raw_files), 0,
                         "raw_pvalues should be deleted without --verbose")

    # -- helpers ----------------------------------------------------------

    @staticmethod
    def _detect_join_column():
        return JOIN_COLUMN


# ===================================================================

if __name__ == "__main__":
    unittest.main()
