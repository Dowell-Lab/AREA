#!/bin/bash 
#SBATCH --job-name=area_run # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=64
#SBATCH --partition long
#SBATCH --mem=500gb # Memory limit
#SBATCH --time=100:00:00 # Time limit hrs:min:sec
#SBATCH --time=90:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/eofiles/area_run.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/eofiles/area_run.%j.err # Standard error log


dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"


#turn on the virtual machine you are using if you are using one#
path_to_venv=$HOME/VENVs/areavenv/
source $path_to_venv/bin/activate

#set paths to AREA and to files to load in
path_to_area=$HOME/AREA/src/
indir=/Shares/down/public/INLCUDE_2024/kallisto_20241030/selfannoated/
commoncolumn=Participant
rank_file=${indir}kallisto_200401lines_participants_normcounts.csv
boolean_attribute_file=${indir}full_HP_binary_attribute.csv
outdirname=$HOME/area_runs/AREA_2025/outdir/
include_sample_file=${indir}include_participants_with_RNA_and_completeT21.csv
include_rank_file_columns=${indir}include_rank_cols_minexp_1.csv
include_boolean_file_columns=${indir}include_bool_cols_min_5_cT21.csv
outdirname_pre=${outdirname}T21_minexp1_mincomobid5T21_HP


echo $rank_file
echo $boolean_attribute_file
echo $outdirname

echo $outdirname_pre

python3 ${path_to_area}AREA_core.py --verbose -od $outdirname_pre -cc $commoncolumn -rf $rank_file -baf $boolean_attribute_file --processes 64 --include_rank_file_columns $include_rank_file_columns --include_boolean_file_columns $include_boolean_file_columns --include_sample_file $include_sample_file

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

