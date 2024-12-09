# AREA

Attribute rank enrichment analysis

## Summary

The goal of AREA is to link binary attributes to value attributes via samples. In this way, we may be able to find value attributes that cause the binary attributes.

In our example, we analyze gene expression data from 254 individuals with Down syndrome, along with their associated medical conditions. Our objective is to identify genes that, when highly expressed, may influence the likelihood of specific medical conditions.  However, it's important to note that, similar to correlation analysis, a connection between a gene and a condition does not imply causation. 

![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

NES stands for Normalized Enrichment Score. A negative NES indicates that high levels of the value (gene) are associated with the true condition of the binary attribute. Conversely, a positive NES suggests that low levels of the value (gene) are associated with the true condition of the binary attribute.

## AREA inputs and outputs

AREA is a tool that takes two CSV files as input. Both CSVs have the same samples, but one CSV has a binary attribute for each sample and one has a value attribute for each sample. The output is a CSV with statistically significant linkages between the binary attributes and the value columns. 

### Input files

Both input CSV files must have a common sample name column. 
In our case, that common sample name column is "Patient".

### Binary Attributes file

One of the CSV files needs to contain binary attributes. In our case, the binary attributes are the disease or disorder associated with each patient, which we call a comorbditity. 
![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/binary_attributes_df.png "binary attributes csv")

### Value file

The other files must have columns that can be ranked by the values within the column. In our example case, the value file has genes and the expression level of those genes in each patient. 

![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/value_df.png "Value csv")

### Output file

The output file has each binary attribute and value type listed in pairs. The rest of the row has the scores for this pair. For instance, at the top of the table below, X and Y are shown. They have a NES that says -Y which means that this gene is higher in the people with Y.
![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

### Filtering the output

The items in the output file have NOT been filtered for significance. To filter for significance you must pick an adjusted p-value column (we include 4 to choose from) and filter out values less than your custom cutoff. 

From the outset, it is unclear which comorbiditys will be linked to which genes in our data set. But, it is likley some genes won't be linked becuase there gene expression is not replicable. Similarly, if a comorbidity is present in all patients or in zero patients, is it unlikely to be significant in AREA. Therefore, we chose not to run these genes and comrbidities thorugh the program. If we did, we would would have to multiple hypothesis correct for each test. ** is this right ** 

## Running AREA

To run AREA, look at the examples in the examples folder. 

The example file slurm_runpsea_fullinputfiles.sh will run AREA on a supercomputer with slurm installed. It will run using the data in test data. You can modify this script to run your own data. 

Often a user will not want to run all data through AREA. As mentioned before, if a binary attribute is true of all samples or no samples it is silly to run AREA. Therefore we have option --include flags. Those flags take in a file like the file "include_binary_attribute_short.csv" which has a list of binary_attribute columns from the binary_attributes file. The rest of the binary_attribute file will be marked as "exclude" in the output file. 
To run the filtered version of AREA use the example file slurm_runpsea_fullinputfiles_filterviainclude.sh.

Finally, we include a run file for simulated data. See more in the description of Simulated data below. 
slurm_runpsea_fullsimulatedfiles.sh

## Simulated data

To avoid genes and comorbidities that cannot be called as significant even after these filtering methods, we decided to create simulated data and discover what patterns AREA is best at finding. 
We created simulated data similar to our data in the notebook_examples//simulateddata-bothdirs.ipynb. This code creates a certain number of simulated genes based on the expression of real genes. We then create comorbiditis that are bais in a sub-group of people with higher or lower expression. 

## Special thank you to the Human Trisomy Project and the INCLUDE Data Coordinating Center

The data used in our example code can be sourced from the INCLUDE Data Hub at https://portal.includedcc.org/, and anyone who would like to use the original data can register for an account on the data hub. Data used from the Human Trisomy Project was converted from kallisto to deseq2 normalized counts using tximport.
