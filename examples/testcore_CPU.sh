#!/bin/bash 
#SBATCH --job-name=area # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=6
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/e_and_o/slurm_test.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/e_and_o/slurm_test.%j.err # Standard error log


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


path_to_venv=$HOME
source $path_to_venv/psea_venv/bin/activate


#set the paths to the files to load in
path_to_area=$HOME/AREA/
indir=$HOME/ChrisO/PSEA/AREA_fast/output/
rank_file=${indir}value_expression.csv
boolean_attribute_file=${indir}comorbid_file.csv
outdirname=$HOME/outarea/initialareaCPU

mkdir -p "$(dirname "$outdirname")"

echo $rank_file
echo $boolean_attribute_file
echo $outdirname

python3 ${path_to_area}run_area.py \
  -od $outdirname \
  -jc Patient \
  -rf $rank_file \
  -bf $boolean_attribute_file \
  -t 60

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"
