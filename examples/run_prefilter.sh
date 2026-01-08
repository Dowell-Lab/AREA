#!/bin/bash 
#SBATCH --job-name=filter # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ozeroff@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/ozeroff/Ozeroff_scratch/ChrisO/AREA_git/chris_AREA/eando/slurm_test.%j.out # Standard output
#SBATCH --error=/Users/ozeroff/Ozeroff_scratch/ChrisO/AREA_git/chris_AREA/eando/slurm_test.%j.err # Standard error log

#turn on the virtual machine you are using if you are using one
path_to_venv=$HOME
source $path_to_venv/jhub_venv/bin/activate 

indir=/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/
outdir=$HOME/Ozeroff_scratch/ChrisO/AREA_git/chris_AREA/output/filtered_data/all_chroms/just_T21/
common_column_name=Participant
values_file=${indir}kallisto_200401lines_participants_normcounts.csv
binary_attribute_file=${indir}full_MONDO_binary_attribute.csv
prefilter_func=$HOME/Ozeroff_scratch/ChrisO/AREA_git/chris_AREA/src/

#chr21 gene file path (needed for running only chr21 arguments, otherwise will n0ot do anything)
chr21_path=$HOME/down_public/INLCUDE_2024/kallisto_20241030/selfannoated/gene_id_no_version_chr21only.csv

# filtered dataframe outputs for AREA
filtered_values_file=${outdir}filtered_values_dataframe.csv
filtered_binary_attribute_file=${outdir}filtered_binary_attributes_dataframe.csv

# reference output lists, with included and excluded genes/comorbids
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

python3 ${prefilter_func}run_prefilter.py \
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
    --patient_comorbid_threshold 0 \
    --min_comorbids_percent 0.05 \
    --max_comorbids_percent 0.95 \
    --min_mean_expression 1.0 \
    --individual_expression_threshold 10 \
    --verbose
##    --t21_only \
##    --t21_column MONDO_complete_trisomy_21 \
## Use these arguments to only run AREA with the complete_trisomy_21 individuals
## Other arguments (Chr21 only, removing cobnfounding comorbids) found in "run_prefilter_py"

# standard output
if [ $? -eq 0 ]; then
    echo "Filtering completed successfully at: $(date '+%Y-%m-%d %H:%M:%S')"
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
