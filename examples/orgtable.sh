indir=$HOME/psea/testdata/
sample_name=Patient
values_file=${indir}value_expression.csv
bianary_attribute_file=${indir}comorbid_file.csv
outdirname=$HOME/outpsea/inputfiles/

mkdir $outdirname


echo $values_file
echo $bianary_attribute_file
echo $outdirname


python3 ../src/psea_organization.py -od $outdirname -sn $sample_name -vf $values_file -baf $bianary_attribute_file




