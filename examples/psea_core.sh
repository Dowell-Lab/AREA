indir=$HOME/outpsea/inputfiles/
#infile=patent_ductus_arteriosus_ENSG00000198743.csv
#infile=eye_disorder_ENSG00000159216.csv
infile=complete_trisomy_21_ENSG00000198743.csv
#infile=eye_disorder_ENSG00000198743.csv
#infile=eye_disorder_ENSG00000201984.csv
#indir=$HOME/psea/testdata/
#infile=testdata.csv
outdir=$HOME/outpsea/
outfilename=$infile

mkdir $outdirname

python3 ../src/psea_core.py -i ${indir}$infile -o ${outdir}$outfilename.pseascore.csv

echo ${outdir}$outfilename.pseascore.csv
