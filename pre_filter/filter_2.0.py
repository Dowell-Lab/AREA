import pandas as pd
import os

def validate_input_files(values_file, binary_attribute_file, sample_name):
    """This function validates if the input files exist and have the required columns."""
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