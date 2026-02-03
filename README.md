# AREA

Attribute rank enrichment algorithum

## Summary

The goal of AREA is to link boolean attributes to rankable attributes via samples. In this way, we may be able to find rankable attributes that cause the boolean attributes or vise versa. Similar to correlation anaysis, causality is not clear from this anaysis alone. Instead these linkages should be considered potientally causal, and downstream experiments should be use to determine causality. 

In our example, we analyze gene expression data from individuals. These individuals have Down syndrome, and the boolean attributes are their associated medical conditions. Our objective is to identify genes that, when highly expressed, may influence the likelihood of specific medical conditions.   

![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

NES stands for Normalized Enrichment Score. A negative NES indicates that low levels of the rankable attribute (in our case gene expression) are associated with the true condition of the boolean attribute. Conversely, a positive NES suggests that high levels of the rankable attribute (in our case gene expression) are associated with the true condition of the boolean attribute.

## AREA inputs and outputs

AREA is a tool that takes two CSV files as input. Both CSVs have the same samples (in our example, individuals), but one CSV has a boolean attributes for each sample and one has a rankable attributes for each sample. The output is a CSV with the statistically significance of the linkages between the boolean columns and the rank columns. 

### Input files

Both input CSV files must have a common sample column. 
In our case, that common sample name column is "Participant".

### Binary Attribute file

One of the CSV files needs to contain boolean attributes. In our case, the boolean attributes are the disease or disorder associated with each patient, which we call a comorbditity. 
![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/boolean_attributes_df.png "boolean attributes csv")

### Value file

The other file must have columns that can be ranked by the ranks within the column. In our example case, the rank file has genes and the expression level of those genes in each patient. 

![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/rank_df.png "Value csv")

### Running pre-filter

In the AREA/pre_filter repository, we have included a filtering bash script called run_prefilter.sh. This script takes in the boolean attribute file and rank file and outputs filtered lists of genes and comorbidities. You can then run AREA either using the filterted subset of genes and comorbidities, or using the full orginal CVS but with the include flags (--include_rank_file_columns, --include_boolean_file_columns, --include_sample_file) which allow you to limit which columns from each of the CSVs are used by AREA.  

In addition to setting designating your input csv files (--boolean_file,--rank_file), you will need to set your out-directory, as well as the common_column_name (name of the column that that is in both the boolean_file and the rank_file).

, the name of the rank file and boolean attribute file, and the names of the files that will be created for included and excluded genes and comorbidities. patient_comorbid_threshold is set to 1, and is the minimum number of comorbidities a patient must have to be included (ensures every patient has at least one comorbid). min_comorbids_percent and max_comorbids_percent are the min and max percent prevalence of a comorbid for the comorbid to be included, and these are set to 0.05 and 0.95 respectively. min_mean_expression is the minimum level of expression a gene must have to be included and is set to 0.1 (in the HTP data example this would be in TPM). individual_expression_threshold is the minimum overall gene expression level a patient must have for the patient to be included in the analysis, and this is set to 10. 

After filtering, the output directory will contain files that have both the included and excluded genes and comorbidities, which can now be used in AREA. A description of each flag can be found in the run_filter.py script which uses argument parser. The individual filtering functions can be found in the filter_functions.py script.

### Output file

The output file has each boolean attribute and rank type listed in pairs. The rest of the row has the scores for this pair. For instance, at the top of the table below, X and Y are shown. They have a NES that says -Y which means that this gene is higher in the people with Y.
![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

### Filtering the output

The items in the output file have NOT been filtered for significance. To filter for significance you must pick an adjusted p-rank column (we include 4 to choose from) and filter out ranks less than your custom cutoff. 

From the outset, it is unclear which comorbiditys will be linked to which genes in our data set. But, it is likley some genes won't be linked becuase there gene expression is not replicable. Similarly, if a comorbidity is present in all patients or in zero patients, is it unlikely to be significant in AREA. Therefore, we chose not to run these genes and comrbidities thorugh the program. If we did, we would would have to multiple hypothesis correct for each test. ** need to change this ** 

## Running AREA

To run AREA, look at the examples in the examples folder. 

The example file slurm_runpsea_fullinputfiles.sh will run AREA on a supercomputer with slurm installed. It will run using the data in test data. You can modify this script to run your own data. 

*** I need to change this*** Often a user will not want to run all data through AREA. As mentioned before, if a boolean attribute is true of all samples or no samples it is silly to run AREA. Therefore we have option --include flags. Those flags take in a file like the file "include_boolean_attribute_short.csv" which has a list of boolean_attribute columns from the boolean_attributes file. The rest of the boolean_attribute file will be marked as "exclude" in the output file. 
To run the filtered version of AREA use the example file slurm_runpsea_fullinputfiles_filterviainclude.sh.

Finally, we include a run file for simulated data. See more in the description of Simulated data below. 
slurm_runpsea_fullsimulatedfiles.sh

## Simulated data

To avoid genes and comorbidities that cannot be called as significant even after these filtering methods, we decided to create simulated data and discover what patterns AREA is best at finding. 
We created simulated data similar to our data in the notebook_examples//simulateddata-bothdirs.ipynb. This code creates a certain number of simulated genes based on the expression of real genes. We then create comorbiditis that are bais in a sub-group of people with higher or lower expression. 

## Special thank you to the Human Trisomy Project and the INCLUDE Data Coordinating Center

The data used in our example code can be sourced from the INCLUDE Data Hub at https://portal.includedcc.org/ and anyone who would like to use the original data can register for an account on the data hub. Data used from the Human Trisomy Project was converted from kallisto to deseq2 normalized counts using tximport.
