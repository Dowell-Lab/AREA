indir=$HOME/psea/testdata/
sample_name=Patient
values_file=${indir}value_expression.csv
bianary_attribute_file=${indir}comorbid_file.csv
outdirname=$HOME/outpsea/
include_values_file=${indir}include_values_short.csv
include_binary_attribute_file=${indir}include_binary_attribute_short.csv

mkdir $outdirname


echo $values_file
echo $bianary_attribute_file
echo $outdirname


python3 ../src/psea_wrapper.py -od $outdirname -sn $sample_name -vf $values_file -baf $bianary_attribute_file --include_values_file $include_values_file --include_bianary_attribute_file $include_binary_attribute_file




