## the original filtering function, not very modular but works ##

import pandas as pd
import os

def run_filtering(patient_comorbid_threshold, 
                  min_comorbids_percent, 
                  max_comorbids_percent, 
                  individual_expression_threshold, 
                  min_mean_expression, 
                  values_file, 
                  binary_attribute_file, 
                  sample_name, 
                  include_values_file, 
                  include_binary_attribute_file, 
                  exclude_values_file, 
                  exclude_binary_attribute_file):
    
    """
    
    This function reads in the cormorbidity designation and gene expression data files, 
    and filters them based on user defined inputs. It then outputs lists of filtered Genes and Comorbidities that can be used for AREA. 
    It requires the following arguments: 
    
    - patient_comorbid_threshold = Threshold for minimum number of comorbidities a patient should have to be included. 

    - min_comorbids_percent = Threshold for the minimum percent prevolence of a comorbidity to retain that comorbidity. 
    
    - max_comorbids_percent = Threshold for the maximum percent prevolence of a comorbidity to retain that comorbidity.
    
    - individual_expression_threshold = the lowest expression value an individual has 
    across all their genes before that individual is thrown out.
    
    - min_mean_expression = Threshold for minimum gene expression level to retain a gene.
    
    - values_file = the file with patient IDs, gene names, and gene expression values.
    
    - binary_attribute_file = the file with patient IDs and comorbidity designations.
    
    - sample_name = the title of the column that has the patient IDs or sample names in the comorbidity and gene expression files.
    
    - include_values_file = Name/location of file to contain filtered list of gene names.
    
    - include_binary_attribute_file = Name/location of file to contain filtered list of cormorbidities.
    
    - exclude_values_file = Name/location of file to contain list of gene names that were excluded based on filtering.
    
    - exclude_binary_attribute_file = Name/location of file to contain list of cormorbidities that were excluded based on filtering.
    
    """
    ###############################
    ## comorbidity csv filtering ##
    ###############################
    
    try:
        # get working dir and read in the CSV file
        file_path = binary_attribute_file
        if not os.path.isfile(binary_attribute_file):
            raise FileNotFoundError(f"Binary attribute file '{binary_attribute_file}' not found.")
        df = pd.read_csv(file_path, index_col=0)
        # check that the sample_name column actually exists
        if sample_name not in df.columns:
            raise ValueError(f"Sample column labeled '{sample_name}' not found in the binary attribute file.")

        # get the number samples
        total_samples = df.shape[0]

        # Filter out patients without comorbidities
        if patient_comorbid_threshold < 0:
            raise ValueError("Patient comorbid threshold must be a positive number.")
        filtered_df = df[df.drop(columns=[sample_name]).sum(axis=1) >= patient_comorbid_threshold]

        # drop sample_names column
        df_noname = filtered_df.drop(columns=[sample_name])

        # Calculate percent occurence of comorbidities
        n_comorbid = df_noname.sum(axis=0).to_frame()
        n_comorbid.columns =["n"]
        n_comorbid["binary_attribute"] = n_comorbid.index
        n_comorbid["percent"] =  n_comorbid["n"]/total_samples

        # Filter based on user input min and max comorbid occurence
        if not (0 <= min_comorbids_percent <= 1 and 0 <= max_comorbids_percent <= 1):
            raise ValueError("Comorbid percentage thresholds must be between 0 and 1.")
        n_comorbid_include = n_comorbid[n_comorbid["percent"]<max_comorbids_percent]
        n_comorbid_include = n_comorbid_include[n_comorbid_include["percent"]>min_comorbids_percent]
        print("Passing comorbidities based on filtering thresholds:")
        print(n_comorbid_include)
        # grab the comorbids that where excluded
        n_comorbid_exclude = n_comorbid[~n_comorbid.index.isin(n_comorbid_include.index)]
        print("Excluded comorbidities based on filtering thresholds:")
        print(n_comorbid_exclude)
        
        # name the output file and save it based on include_binary_attribute_file argument input
        if not os.access(os.path.dirname(include_binary_attribute_file), os.W_OK):
            raise PermissionError(f"Cannot write to directory '{include_binary_attribute_file}'. Did you designate the ouput directoty and filename correctly?")
        outdir = include_binary_attribute_file
        n_comorbid_include[["binary_attribute"]].to_csv(outdir, header=False, index=False)
        print("INCLUDED COMORBIDS FILE SAVED")
        
        # do the same for the exclude file
        if not os.access(os.path.dirname(exclude_binary_attribute_file), os.W_OK):
            raise PermissionError(f"Cannot write to directory '{exclude_binary_attribute_file}'. Did you designate the ouput directoty and filename correctly for the exclude file?")
        outdir = exclude_binary_attribute_file
        n_comorbid_exclude[["binary_attribute"]].to_csv(outdir, header=False, index=False)
        print("EXCLUDED COMORBIDS FILE ALSO SAVED")
        
    except Exception as e:
        print(f"Error occurred during comorbidity filtering: {e}")
    
    ###################################
    ## gene expression csv filtering ##
    ###################################
    
    try:
        # get working dir and read in the CSV file for gene expression data
        file_path = values_file
        if not os.path.isfile(values_file):
            raise FileNotFoundError(f"Values file '{values_file}' not found.")
        df = pd.read_csv(file_path, index_col=0)
        
        # Validate if the sample_name column exists
        if sample_name not in df.columns:
            raise ValueError(f"Sample column '{sample_name}' not found in the values file.")
        # drop samples column
        expression_df_nosamples = df.drop(columns=[sample_name])

        # Remove individuals (rows) where all gene expression values are below the threshold
        individual_mask = expression_df_nosamples.apply(lambda row: (row >= individual_expression_threshold).any(), axis=1)
        expression_df_nosamples = expression_df_nosamples[individual_mask]

        # get mean expression, turn into pandas dataframe
        meansdf=expression_df_nosamples.mean(axis=0).to_frame()
        meansdf.columns = ["mean_value"]
        meansdf["valuename"]= meansdf.index
        
        # Validate if min_mean_expression is a valid threshold
        if min_mean_expression < 0:
            raise ValueError("Minimum mean expression threshold must be non-negative.")

        # filter based on mean expression
        meansdf_include = meansdf[meansdf["mean_value"]>min_mean_expression]
        print("Passing genes based on mean threshold for expression:")
        print(meansdf_include)
        # get the excluded genes
        meansdf_exclude = meansdf[~meansdf.index.isin(meansdf_include.index)]
        print("Excluded genes based on mean threshold for expression:")
        print(meansdf_exclude)

        # save output file based on include_values_file argument input
        if not os.access(os.path.dirname(include_values_file), os.W_OK):
            raise PermissionError(f"Cannot write to directory '{include_values_file}'. Did you designate the ouput directoty and filename correctly?")
       
        outdir = include_values_file
        meansdf_include[["valuename"]].to_csv(outdir, header=False, index=False)
        print("INCLUDE EXPRESSION FILE SAVED")
        
        # do the same for the excluded genes
        if not os.access(os.path.dirname(exclude_values_file), os.W_OK):
            raise PermissionError(f"Cannot write to directory '{exclude_values_file}'. Did you designate the ouput directoty and filename correctly for the excluded genes file?")
       
        outdir = exclude_values_file
        meansdf_exclude[["valuename"]].to_csv(outdir, header=False, index=False)
        print("EXCLUDED EXPRESSION FILE ALSO SAVED")
    
    except Exception as e:
        print(f"Error occurred during gene expression filtering: {e}")
        
    ####Some ideas going forward:
    # Need to test all of the error handling, but one so far look good
    