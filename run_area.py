#!/usr/bin/env python
"""
Top-level entry point for AREA.

Usage:
    python run_area.py --config area_config.yaml
    python run_area.py -bf bools.csv -rf ranks.csv -jc sample_id -od results/
"""

from area.cli import main

if __name__ == "__main__":
    main()
