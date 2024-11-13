#!/bin/bash 
#SBATCH --job-name=filter_psea # Job name
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

indir=../testdata/
outdir=$HOME/ChrisO/PSEA_OUTPUT/
sample_name=Patient
values_file=${indir}value_expression.csv
binary_attribute_file=${indir}comorbid_file.csv
include_values_file=${outdir}include_values_long.csv
include_binary_attribute_file=${outdir}include_binary_attribute_long.csv
exclude_values_file=${outdir}exclude_values_long.csv
exclude_binary_attribute_file=${outdir}exclude_binary_attribute_long.csv

mkdir -p $outdir

echo "Starting filtering at: $(date '+%Y-%m-%d %H:%M:%S') ..."
echo "Filtering the following files:"
echo $values_file
echo $binary_attribute_file
echo "Output files:"
echo $include_values_file
echo $include_binary_attribute_file
echo "Filtering running..."

python3 run_filter.py \
    -sn $sample_name \
    -vf $values_file \
    -baf $binary_attribute_file \
    --include_values_file $include_values_file \
    --exclude_values_file $exclude_values_file \
    --include_binary_attribute_file $include_binary_attribute_file \
    --exclude_binary_attribute_file $exclude_binary_attribute_file \
    --patient_comorbid_threshold 1 \
    --min_comorbids_percent 0.05 \
    --max_comorbids_percent 0.95 \
    --min_mean_expression 0.1 \
    --individual_expression_threshold 10


# In standard output check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "Filtering completed successfully at: $(date '+%Y-%m-%d %H:%M:%S')"
else
    echo "Filtering failed at: $(date '+%Y-%m-%d %H:%M:%S')"
fi