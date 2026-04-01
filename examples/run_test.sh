#!/bin/bash 
#SBATCH --job-name=area_run # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --partition long
#SBATCH --mem=500gb # Memory limit
#SBATCH --time=100:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/eofiles/area_run.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/eofiles/area_run.%j.err # Standard error log


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#turn on the virtual machine you are using if you are using one#
path_to_venv=$HOME/VENVs/areavenv39/
source $path_to_venv/bin/activate

path_to_area=$HOME/AREA/

python ${path_to_area}run_tests.py -v

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"
