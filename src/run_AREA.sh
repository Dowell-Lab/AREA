#!/bin/bash 
#SBATCH --job-name=area_filtered # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ozeroff@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=190
#SBATCH --partition highmem
#SBATCH --mem=190gb # Memory limit
#SBATCH --time=100:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/ozeroff/ChrisO/eando/slurm_test.%j.out # Standard output
#SBATCH --error=/Users/ozeroff/ChrisO/eando/slurm_test.%j.err # Standard error log


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#turn on the virtual machine you are using if you are using one#
path_to_venv=$HOME
source $path_to_venv/jhub_venv/bin/activate


#set paths to AREA and to files to load in
path_to_area=$HOME/ChrisO/PSEA/AREA_fast/src/
indir=$HOME/ChrisO/PSEA/AREA_fast/output/filtered_gene_comorbids/all_chrom/
commoncolumn=Patient
values_file=${indir}filtered_values_dataframe.csv
bianary_attribute_file=${indir}filtered_binary_attributes_dataframe.csv
outdirname=$HOME/ChrisO/PSEA/AREA_fast/output/AREA_output/all_chrom/

mkdir -p "$outdirname"

echo $values_file
echo $bianary_attribute_file
echo $outdirname


python3 ${path_to_area}AREA_core.py -od $outdirname -cc $commoncolumn -vf $values_file -baf $bianary_attribute_file --processes 190

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"
