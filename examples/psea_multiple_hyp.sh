indir=$HOME/outpsea/
infilename=psea_scores_20240918-172357.csv
outdir=$HOME/outpsea/
outfilename=${infilename}.adjpval.csv

source /Users/allenma/psea_venv/bin/activate

mkdir $outdirname

python3 ../src/psea_apval.py -i ${indir}$infilename -o ${outdir}$outfilename

echo ${outdir}$outfilename 
