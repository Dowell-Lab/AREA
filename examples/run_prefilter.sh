#!/bin/bash 
#SBATCH --job-name=filter # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ozeroff@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/temp/ChrisO/eando/slurm_test.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/temp/ChrisO/eando/slurm_test.%j.err # Standard error log

#turn on the virtual machine you are using if you are using one
path_to_venv=$HOME
source $path_to_venv/jhub_venv/bin/activate 

indir=~/ChrisO/PSEA_OUTPUT/
outdir=$HOME/ChrisO/PSEA/AREA_fast/output/filtered_gene_comorbids/all_chrom/
common_column_name=Patient
values_file=${indir}first400_express_values_df.csv
binary_attribute_file=${indir}first400_comorbid_df.csv

### Main outputs - filtered dataframes for AREA
filtered_values_file=${outdir}filtered_values_dataframe.csv
filtered_binary_attribute_file=${outdir}filtered_binary_attributes_dataframe.csv

### Reference output lists, with included and excluded genes/comorbids
include_values_file=${outdir}include_values_long.csv
include_binary_attribute_file=${outdir}include_binary_attribute_long.csv
exclude_values_file=${outdir}exclude_values_long.csv
exclude_binary_attribute_file=${outdir}exclude_binary_attribute_long.csv

mkdir -p $outdir

echo "Starting filtering at: $(date '+%Y-%m-%d %H:%M:%S') ..."
echo "Filtering the following files:"
echo $values_file
echo $binary_attribute_file
echo ""
echo "Main output files (filtered dataframes for AREA):"
echo $filtered_values_file
echo $filtered_binary_attribute_file
echo ""
echo "Reference output files (lists):"
echo $include_values_file
echo $include_binary_attribute_file
echo ""
echo "Filtering running..."

python3 run_prefilter.py \
    -cc $common_column_name \
    -vf $values_file \
    -baf $binary_attribute_file \
    -od $outdir \
    --filtered_values_file $filtered_values_file \
    --filtered_binary_attribute_file $filtered_binary_attribute_file \
    --include_values_file $include_values_file \
    --exclude_values_file $exclude_values_file \
    --include_binary_attribute_file $include_binary_attribute_file \
    --exclude_binary_attribute_file $exclude_binary_attribute_file \
    --patient_comorbid_threshold 1 \
    --min_comorbids_percent 0.05 \
    --max_comorbids_percent 0.95 \
    --min_mean_expression 1.0 \
    --individual_expression_threshold 10 \
    --verbose
    
# 0.1 was the original for min mean, trying 1.0
# check if the script completed successfully - standard output
if [ $? -eq 0 ]; then
    echo "Filtering completed successfully at: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
    echo "=== NEXT STEPS ==="
    echo "Use these filtered dataframes with AREA:"
    echo "python3 area_core.py \\"
    echo "  -baf $filtered_binary_attribute_file \\"
    echo "  -vf $filtered_values_file \\"
    echo "  -cc $common_column_name \\"
    echo "  -od your_area_output_dir/ \\"
    echo "  --processes 60"
else
    echo "Filtering failed at: $(date '+%Y-%m-%d %H:%M:%S')"
fi