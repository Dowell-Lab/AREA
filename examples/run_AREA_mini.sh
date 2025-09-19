#!/bin/bash 
#SBATCH --job-name=area_run # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --partition long
#SBATCH --mem=500gb # Memory limit
#SBATCH --time=100:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/eofiles/area_test.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/eofiles/area_test.%j.err # Standard error log


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#turn on the virtual machine you are using if you are using one#
path_to_venv=$HOME/VENVs/areavenv/
source $path_to_venv/bin/activate


#set paths to AREA and to files to load in
path_to_area_main=$HOME/AREA/
path_to_area_testfiles=${path_to_area_main}testdata/
path_to_area=${path_to_area_main}src/
indir=/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/
commoncolumn=Participant
rank_file=${indir}miniRNAvalues_participants_normcounts.csv
boolean_attribute_file=${indir}mini_HP_binary_attribute.csv
outdirname=$HOME/area_runs/AREA_2025/outdir/mini

mkdir -p "$outdirname"

echo $rank_file
echo $boolean_attribute_file
echo $outdirname

echo $outdirname

python3 ${path_to_area}AREA_core.py --verbose -od $outdirname -cc $commoncolumn -rf $rank_file -baf $boolean_attribute_file --processes 64

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

