#!/bin/bash 
#SBATCH --job-name=psea # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/e_and_o/slurm_test.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/e_and_o/slurm_test.%j.err # Standard error log

#turn on the virtual machine you are using if you are using one
path_to_venv=$HOME
source $path_to_venv/psea_venv/bin/activate

#set the paths to the files to load in
path_to_psea=$HOME
indir=$HOME/psea/testdata/
sample_name=Patient
values_file=${indir}value_expression.csv
bianary_attribute_file=${indir}comorbid_file.csv
outdirname=$HOME/outpsea/

mkdir $outdirname


echo $values_file
echo $bianary_attribute_file
echo $outdirname


python3 ${path_to_psea}/psea/src/psea_wrapper.py -od $outdirname -sn $sample_name -vf $values_file -baf $bianary_attribute_file --processes 60




