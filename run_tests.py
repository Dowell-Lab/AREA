#!/usr/bin/env python
"""
Run the AREA test suite.

Usage:
    python run_tests.py
    python run_tests.py -v       # verbose output
"""

import sys
import unittest


if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir="tests", pattern="test_*.py")

    verbosity = 2 if "-v" in sys.argv or "--verbose" in sys.argv else 1
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)

    sys.exit(0 if result.wasSuccessful() else 1)
