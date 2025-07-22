#!/usr/bin/env python3
"""
Prefilter script for AREA analysis
Outputs filtered dataframes for use with AREA
"""

import argparse
import os
import time
from filter_functions import run_filtering

def main():
    """Main function for standalone prefiltering"""
    parser = argparse.ArgumentParser(description="Prefilter data for AREA analysis")

    # Required arguments
    parser.add_argument('-vf', '--values_file', type=str, required=True,
                        help='File with patient IDs, gene names, and expression values')
    parser.add_argument('-baf', '--binary_attribute_file', type=str, required=True,
                        help='File with patient IDs and comorbidity designations')
    parser.add_argument('-cc', '--common_column_name', type=str, required=True,
                        help='Column name that contains patient IDs/names')
    parser.add_argument('-od', '--outdir', type=str, required=True,
                        help='Output directory for filtered files')

    # Filtering thresholds
    parser.add_argument('--patient_comorbid_threshold', type=int, default=1,
                        help='Minimum number of comorbidities a patient should have (default: 1)')
    parser.add_argument('--min_comorbids_percent', type=float, default=0.05,
                        help='Minimum percent prevalence of comorbidity to retain (default: 0.05)')
    parser.add_argument('--max_comorbids_percent', type=float, default=0.95,
                        help='Maximum percent prevalence of comorbidity to retain (default: 0.95)')
    parser.add_argument('--min_mean_expression', type=float, default=1.0,
                        help='Minimum mean gene expression to retain gene (default: 1.0)')
    parser.add_argument('--individual_expression_threshold', type=int, default=10,
                        help='Minimum individual expression threshold (default: 10)')

    # Output file specifications - UPDATED FOR DATAFRAMES
    parser.add_argument('--filtered_values_file', type=str, default=None,
                        help='Output file for filtered gene expression dataframe (default: auto-generated)')
    parser.add_argument('--filtered_binary_attribute_file', type=str, default=None,
                        help='Output file for filtered comorbidity dataframe (default: auto-generated)')
    
    # Optional: still output the old format lists for reference
    parser.add_argument('--include_values_file', type=str, default=None,
                        help='Output file for list of genes to include (default: auto-generated)')
    parser.add_argument('--include_binary_attribute_file', type=str, default=None,
                        help='Output file for list of comorbidities to include (default: auto-generated)')
    parser.add_argument('--exclude_values_file', type=str, default=None,
                        help='Output file for list of genes to exclude (default: auto-generated)')
    parser.add_argument('--exclude_binary_attribute_file', type=str, default=None,
                        help='Output file for list of comorbidities to exclude (default: auto-generated)')

    # Options
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')

    args = parser.parse_args()

    # Validate input files
    if not os.path.isfile(args.values_file):
        raise FileNotFoundError(f"Values file not found: {args.values_file}")
    if not os.path.isfile(args.binary_attribute_file):
        raise FileNotFoundError(f"Binary attribute file not found: {args.binary_attribute_file}")

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Generate timestamp for output files if not specified
    timestr = time.strftime("%Y%m%d-%H%M%S")

    # Set default output file names for FILTERED DATAFRAMES (main outputs)
    if args.filtered_values_file is None:
        args.filtered_values_file = os.path.join(args.outdir, f"filtered_values_{timestr}.csv")
    if args.filtered_binary_attribute_file is None:
        args.filtered_binary_attribute_file = os.path.join(args.outdir, f"filtered_binary_attributes_{timestr}.csv")

    # Set default output file names for LISTS (reference outputs)
    if args.include_values_file is None:
        args.include_values_file = os.path.join(args.outdir, f"include_values_list_{timestr}.csv")
    if args.include_binary_attribute_file is None:
        args.include_binary_attribute_file = os.path.join(args.outdir, f"include_binary_attributes_list_{timestr}.csv")
    if args.exclude_values_file is None:
        args.exclude_values_file = os.path.join(args.outdir, f"exclude_values_list_{timestr}.csv")
    if args.exclude_binary_attribute_file is None:
        args.exclude_binary_attribute_file = os.path.join(args.outdir, f"exclude_binary_attributes_list_{timestr}.csv")

    # Print configuration if verbose
    if args.verbose:
        print("=== PREFILTER CONFIGURATION ===")
        print(f"Input files:")
        print(f"  Values: {args.values_file}")
        print(f"  Binary attributes: {args.binary_attribute_file}")
        print(f"  Common column: {args.common_column_name}")
        print(f"")
        print(f"Filtering thresholds:")
        print(f"  Patient comorbid threshold: {args.patient_comorbid_threshold}")
        print(f"  Comorbidity prevalence range: {args.min_comorbids_percent} - {args.max_comorbids_percent}")
        print(f"  Minimum mean expression: {args.min_mean_expression}")
        print(f"  Individual expression threshold: {args.individual_expression_threshold}")
        print(f"")
        print(f"Main output files (filtered dataframes):")
        print(f"  Filtered values dataframe: {args.filtered_values_file}")
        print(f"  Filtered binary attributes dataframe: {args.filtered_binary_attribute_file}")
        print(f"")
        print(f"Reference output files (lists):")
        print(f"  Include values list: {args.include_values_file}")
        print(f"  Include binary attributes list: {args.include_binary_attribute_file}")
        print(f"  Exclude values list: {args.exclude_values_file}")
        print(f"  Exclude binary attributes list: {args.exclude_binary_attribute_file}")
        print(f"")

    print(f"Starting filtering at: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    try:
        # Run the filtering
        run_filtering(
            patient_comorbid_threshold=args.patient_comorbid_threshold,
            min_comorbids_percent=args.min_comorbids_percent,
            max_comorbids_percent=args.max_comorbids_percent,
            individual_expression_threshold=args.individual_expression_threshold,
            min_mean_expression=args.min_mean_expression,
            values_file=args.values_file,
            binary_attribute_file=args.binary_attribute_file,
            sample_name=args.common_column_name,
            filtered_values_file=args.filtered_values_file,
            filtered_binary_attribute_file=args.filtered_binary_attribute_file,
            include_values_file=args.include_values_file,
            include_binary_attribute_file=args.include_binary_attribute_file,
            exclude_values_file=args.exclude_values_file,
            exclude_binary_attribute_file=args.exclude_binary_attribute_file
        )

        print(f"Filtering completed successfully at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"")
        print(f"=== MAIN OUTPUTS (Use these with AREA) ===")
        print(f"Filtered gene expression dataframe: {args.filtered_values_file}")
        print(f"Filtered comorbidity dataframe: {args.filtered_binary_attribute_file}")
        print(f"")
        print(f"You can now run AREA analysis using these filtered dataframes:")
        print(f"")
        print(f"python3 area_core.py \\")
        print(f"  -baf {args.filtered_binary_attribute_file} \\")
        print(f"  -vf {args.filtered_values_file} \\")
        print(f"  -cc {args.common_column_name} \\")
        print(f"  -od your_output_dir/ \\")
        print(f"  --processes 60")

    except Exception as e:
        print(f"Filtering failed at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Error: {e}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())