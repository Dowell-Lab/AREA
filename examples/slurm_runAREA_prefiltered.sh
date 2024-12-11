#!/bin/bash 
#SBATCH --job-name=AREA # Job name
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

indir=$HOME/ChrisO/PSEA_OUTPUT/
sample_name=Patient
values_file=${indir}value_expression.csv
bianary_attribute_file=${indir}comorbid_file.csv
outdirname=$HOME/ChrisO/AREA_OUTPUT/
include_values_file=${indir}include_values_long.csv # in this example we are looking at the included genes
include_binary_attribute_file=${indir}include_binary_attribute_long.csv # and binary atribute file values after filtering, not the excluded ones
#exclude_values_file=${indir}exclude_values_long.csv
#exclude_binary_attribute_file=${indir}exclude_binary_attribute_long.csv


mkdir $outdirname

echo $values_file
echo $bianary_attribute_file
echo $outdirname


python3 ../src/psea_wrapper.py \
 -od $outdirname \
 -sn $sample_name \
 -vf $values_file \
 -baf $bianary_attribute_file \
 -ivf $include_values_file \
 -ibaf $include_binary_attribute_file \
 --processes 60


