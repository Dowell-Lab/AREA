# PSEA
Participant set enrichment analysis

## Summary
The goal of PSEA is to link binary attributes to values attributes via samples. In this way we may be able to find value  attributes that cause the binary attributes.

In our example, we analyze gene expression data from 254 individuals with Down syndrome, along with their associated medical conditions. Our objective is to identify genes that, when highly expressed, may influence the likelihood of specific medical conditions.  However, it's important to note that, similar to correlation analysis, a connection between a gene and a condition does not imply causation. 

![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

NES stands for Normalized Enrichment Score. A negative NES indicates that high levels of the value (gene) are associated with the true condition of the binary attribute. Conversely, a positive NES suggests that low levels of the value (gene) are associated with the true condition of the binary attribute.

## PSEA inputs and outputs
PSEA is a tool that takes two csv files as input. Both csvs have the same samples, but one csv has a biarary attirubute for each sample and one has a value attibutes for each sample. The output is a csv with statically signigant linkages between the biarary attirbutes and the value columns. 

### Input files
Both input csv files must have a comman samplename column. 
In our case that common samplename is "Patient".

### Binary Attributes file
One of the csv files needs to have in it bianary attributes. In our case the bianrary attributes are the disease/disorder of the patients called comorbdititys. 
![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/binary_attributes_df.png "binary attributes csv")

### Value file

One of the files must have columns that can be ranked by the values within the column. In the example case the value file has genes and the expression level of those genes in each Patient. 

![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/value_df.png "Value csv")

### Output file
The output file has each binary attribute and value type listed in pairs. The rest of the row has the scores for this pair. For instance, at the top of the table below, X and Y are shown. They have a NES that says -Y which means that this gene is higher in the people with Y.
![alt text](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

### Filtering the output
 The items in the output file have NOT been filtered for significance. To filter for significance you must pick a adjusted p-value coulm (we include 4 to chose from) and fileter on values less then your custom cutoff. 

It is unclear which comorbiditys will be linked to which genes in our data set. But it is likley some genes won't be linked becuase there gene expression is not replicable. Similarly if a comorbidity is true in all people or no people is it unlikely to be significant in PSEA. Therefore we chose not to run those gene/comrbidityes thorugh the program, as if we run them we must multiple hyptoehsis correct for each test. 

## Running PSEA

To run PSEA, look at the examples in the examples folder. 

THe example file slurm_runpsea_fullinputfiles.sh will run PSEA on a super computer with slurm installed. It will run using the data in test data. You can modify this script to run your own data. 

Often a user will not want to run all data through PSEA. For instance if a bianary attribute is true of all samples or no samples it is silly to run PSEA. Therefore we have option --include flags. Those flags take in a file like the file "include_binary_attribute_short.csv" which has a list of binary_attribute columns from the binary_attributes file. The rest of the binary_attribute file will be marked as "exclude" in the output file. 
To run the filtered version of PSEA use the example file slurm_runpsea_fullinputfiles_filterviainclude.sh.

Finnally we also include a run file for simulated data. See more in the description of Simulated data below. 
slurm_runpsea_fullsimulatedfiles.sh

## Simulated data

To aviod genes and comorbdities that can not be called as significant anyway, we decided to create simulated data and discover what pattens PSEA is best at finding. 
We created simulated data similar to our data in the notebook_examples//simulateddata-bothdirs.ipynb. This code creates a certain number of simulated genes based on the expression of real genes. We then create comorbiditeis that are bais in a sub group of people with higher or lower expression. 


