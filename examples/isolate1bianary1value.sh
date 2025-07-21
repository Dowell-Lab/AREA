indir=$HOME/psea/testdata/
sample_name=Patient
values_file=${indir}value_expression.csv
value_name=ENSG00000198743
#value_name=ENSG00000201984
#value_name=ENSG00000159216
bianary_attribute_file=${indir}comorbid_file.csv
#bianary_attribute_name=patent_ductus_arteriosus
#bianary_attribute_name=eye_disorder
bianary_attribute_name=complete_trisomy_21
outdirname=$HOME/outpsea/inputfiles/
outfilename=${outdirname}${bianary_attribute_name}_${value_name}.csv

mkdir $outdirname


echo $values_file
echo $bianary_attribute_file
echo $outfilename


python3 ../src/make_rankable_file.py -of $outfilename -sn $sample_name -vf $values_file -baf $bianary_attribute_file -vn $value_name -ban $bianary_attribute_name

echo $outfilename



