import pandas as pd
import os

def validate_input_files(values_file, binary_attribute_file, sample_name):
    """ This function reads in the input csv files as dataframes, 
    and validates if the input files exist and have the required 'sample_name' columns. """
    
    if not os.path.isfile(binary_attribute_file):
        raise FileNotFoundError(f"Binary attribute file '{binary_attribute_file}' not found.")
    if not os.path.isfile(values_file):
        raise FileNotFoundError(f"Values file '{values_file}' not found.")
    
    binary_df = pd.read_csv(binary_attribute_file, index_col=0)
    values_df = pd.read_csv(values_file, index_col=0)
    
    if sample_name not in binary_df.columns:
        raise ValueError(f"Sample column labeled '{sample_name}' not found in the binary attribute file.")
    if sample_name not in values_df.columns:
        raise ValueError(f"Sample column labeled '{sample_name}' not found in the values file.")
    
    return binary_df, values_df

def filter_comorbidities(binary_df, patient_comorbid_threshold, min_comorbids_percent, max_comorbids_percent, sample_name):
    """ This function filters the comorbidity data based on the user defined thresholds. """
    total_samples = binary_df.shape[0]

    if patient_comorbid_threshold < 0:
        raise ValueError("Patient comorbid threshold must be a positive number.")
    if not (0 <= min_comorbids_percent <= 1 and 0 <= max_comorbids_percent <= 1):
        raise ValueError("Comorbid percentage thresholds must be between 0 and 1.")

    # Filter out patients without comorbidities
    filtered_df = binary_df[binary_df.drop(columns=[sample_name]).sum(axis=1) >= patient_comorbid_threshold]
    df_noname = filtered_df.drop(columns=[sample_name])

    # Calculate percent occurrence of comorbidities
    n_comorbid = df_noname.sum(axis=0).to_frame()
    n_comorbid.columns = ["n"]
    n_comorbid["binary_attribute"] = n_comorbid.index
    n_comorbid["percent"] = n_comorbid["n"]/total_samples

    # Filter comorbidities based on prevalence thresholds
    n_comorbid_include = n_comorbid[n_comorbid["percent"]<max_comorbids_percent]
    n_comorbid_include = n_comorbid_include[n_comorbid_include["percent"]>min_comorbids_percent]
    n_comorbid_exclude = n_comorbid[~n_comorbid.index.isin(n_comorbid_include.index)]
    
    return n_comorbid_include, n_comorbid_exclude

def filter_gene_expression(values_df, individual_expression_threshold, min_mean_expression, sample_name):
    """ This function filters the gene expression data based on the user defined thresholds. """
    expression_df_nosamples = values_df.drop(columns=[sample_name])

    # Filter individuals (rows) where all gene expression values are below the threshold
    individual_mask = expression_df_nosamples.apply(lambda row: (row >= individual_expression_threshold).any(), axis=1)
    expression_df_nosamples = expression_df_nosamples[individual_mask]

    # Calculate mean expression
    meansdf = expression_df_nosamples.mean(axis=0).to_frame()
    meansdf.columns = ["mean_value"]
    meansdf["valuename"] = meansdf.index

    if min_mean_expression < 0:
        raise ValueError("Minimum mean expression threshold must be non-negative.")

    # Filter genes based on mean expression threshold
    meansdf_include = meansdf[meansdf["mean_value"] > min_mean_expression]
    meansdf_exclude = meansdf[~meansdf.index.isin(meansdf_include.index)]

    return meansdf_include, meansdf_exclude

def save_filtered_data(n_comorbid_include, 
                       n_comorbid_exclude, 
                       include_binary_attribute_file, 
                       exclude_binary_attribute_file, 
                       meansdf_include, 
                       meansdf_exclude, 
                       include_values_file, 
                       exclude_values_file):
    """ This function saves the filtered include/exclude data to CSV files in the user defined directory. """
    # Save include comorbid data
    if not os.access(os.path.dirname(include_binary_attribute_file), os.W_OK):
        raise PermissionError(f"Cannot write to directory '{include_binary_attribute_file}'. Did you designate the ouput directoty and filename correctly?.")
    outdir = include_binary_attribute_file
    n_comorbid_include[["binary_attribute"]].to_csv(outdir, header=False, index=False)
    print(f"Included comorbid data saved to {include_binary_attribute_file}")
    # Save excluded comorbid data
    if not os.access(os.path.dirname(exclude_binary_attribute_file), os.W_OK):
        raise PermissionError(f"Cannot write to directory '{exclude_binary_attribute_file}'. Did you designate the ouput directoty and filename correctly?.")
    outdir = exclude_binary_attribute_file
    n_comorbid_exclude[["binary_attribute"]].to_csv(outdir, header=False, index=False)
    print(f"Excluded comorbid data saved to {exclude_binary_attribute_file}")
    
    #save included expression data
    if not os.access(os.path.dirname(include_values_file), os.W_OK):
        raise PermissionError(f"Cannot write to directory '{include_values_file}'. Did you designate the ouput directoty and filename correctly?")
    outdir = include_values_file
    meansdf_include[["valuename"]].to_csv(outdir, header=False, index=False)
    print(f"Included comorbid data saved to {include_values_file}")
    #save excluded expression data
    if not os.access(os.path.dirname(exclude_values_file), os.W_OK):
        raise PermissionError(f"Cannot write to directory '{exclude_values_file}'. Did you designate the ouput directoty and filename correctly for the excluded genes file?")
    outdir = exclude_values_file
    meansdf_exclude[["valuename"]].to_csv(outdir, header=False, index=False)
    print(f"Excluded comorbid data saved to {exclude_values_file}")

    
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
    
    """ This is the main function that runs comorbidity and gene expression filtering. """
    try:
        # Validate input files and load data
        binary_df, values_df = validate_input_files(values_file, 
                                                    binary_attribute_file, 
                                                    sample_name)

        # Comorbidity filtering
        n_comorbid_include, n_comorbid_exclude = filter_comorbidities(binary_df, 
                                                                      patient_comorbid_threshold, 
                                                                      min_comorbids_percent, 
                                                                      max_comorbids_percent, 
                                                                      sample_name)
        
        print("Passing comorbidities based on filtering thresholds:")
        print(n_comorbid_include)
        print("Excluded comorbidities based on filtering thresholds:")
        print(n_comorbid_exclude)

        # Gene expression filtering
        meansdf_include, meansdf_exclude = filter_gene_expression(values_df, 
                                                                  individual_expression_threshold, 
                                                                  min_mean_expression, 
                                                                  sample_name)
        
        print("Passing genes based on mean threshold for expression:")
        print(meansdf_include)
        print("Excluded genes based on mean threshold for expression:")
        print(meansdf_exclude)
        
        # Save filtered comorbidities
        save_filtered_data(n_comorbid_include, 
                           n_comorbid_exclude, 
                           include_binary_attribute_file, 
                           exclude_binary_attribute_file, 
                           meansdf_include, 
                           meansdf_exclude, 
                           include_values_file, 
                           exclude_values_file)
    
    except Exception as e:
        print(f"Error occurred during filtering: {e}")