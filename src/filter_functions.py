import pandas as pd
import os

def load_chr21_genes(chr21_file):
    """Load chromosome 21 gene list and return base Ensembl IDs (without version)"""
    if chr21_file and os.path.isfile(chr21_file):
        # Read file with no header - first line is data
        chr21_df = pd.read_csv(chr21_file, header=None, names=['gene_id'])
        # These IDs don't have versions, so they're already base IDs
        chr21_genes = set(chr21_df['gene_id'].tolist())
        print(f"  Loaded {len(chr21_genes)} chromosome 21 genes")
        return chr21_genes
    return None

def filter_genes_by_chromosome(genes_list, chr21_genes, chr21_only=True):
    """Filter genes based on chromosome 21 membership"""
    if chr21_genes is None:
        return genes_list  # No filtering if no chr21 file provided
    
    filtered = []
    for gene in genes_list:
        # Extract base ID (part before the period)
        base_id = gene.split('.')[0] if '.' in gene else gene
        
        # Check if gene is on chr21
        is_chr21 = base_id in chr21_genes
        
        # Include based on chr21_only flag
        if (chr21_only and is_chr21) or (not chr21_only and not is_chr21):
            filtered.append(gene)
    
    return filtered

def filter_t21_patients(binary_df, values_df, t21_column, sample_name):
    """Filter to keep only patients with complete trisomy 21"""
    if t21_column not in binary_df.columns:
        raise ValueError(f"T21 column '{t21_column}' not found in binary attributes file. "
                        f"Available columns include: {list(binary_df.columns[:10])}...")
    
    # Convert string booleans if needed
    if binary_df[t21_column].dtype == 'object':
        binary_df[t21_column] = binary_df[t21_column].map({'True': True, 'False': False, 1: True, 0: False})
    
    # Get T21 patients based on the Participant column
    t21_participants = binary_df[binary_df[t21_column] == True][sample_name].unique()
    print(f"  Found {len(t21_participants)} patients with {t21_column}")
    
    # Filter both dataframes using the Participant column
    binary_df_filtered = binary_df[binary_df[sample_name].isin(t21_participants)]
    values_df_filtered = values_df[values_df[sample_name].isin(t21_participants)]
    
    print(f"  After T21 filtering:")
    print(f"    Binary attributes: {binary_df_filtered.shape}")
    print(f"    Values: {values_df_filtered.shape}")
    
    # Verify alignment
    assert len(binary_df_filtered[sample_name].unique()) == len(values_df_filtered[sample_name].unique()), \
        "Filtered dataframes must have same number of unique participants!"
    
    return binary_df_filtered, values_df_filtered

def validate_input_files(values_file, binary_attribute_file, sample_name, t21_only=False, t21_column=None):
    """ This function reads in the input csv files as dataframes, 
    and validates if the input files exist and have the required structure. """
    
    if not os.path.isfile(binary_attribute_file):
        raise FileNotFoundError(f"Binary attribute file '{binary_attribute_file}' not found.")
    if not os.path.isfile(values_file):
        raise FileNotFoundError(f"Values file '{values_file}' not found.")
    
    # Read with index_col=0 to use the first column as index
    binary_df = pd.read_csv(binary_attribute_file, index_col=0)
    values_df = pd.read_csv(values_file, index_col=0)
    
    print(f"Initial data shapes:")
    print(f"  Binary attributes: {binary_df.shape}")
    print(f"  Values: {values_df.shape}")
    
    # Check if Participant column exists
    if sample_name not in binary_df.columns or sample_name not in values_df.columns:
        raise ValueError(f"Sample column labeled '{sample_name}' not found in one or both files.")
    
    print(f"Found '{sample_name}' column in both files")
    
    # Find common participants based on the Participant column values
    participants_binary = set(binary_df[sample_name].unique())
    participants_values = set(values_df[sample_name].unique())
    common_participants = participants_binary & participants_values
    
    print(f"  Binary file has {len(participants_binary)} unique participants")
    print(f"  Values file has {len(participants_values)} unique participants")
    print(f"  Found {len(common_participants)} common participants")
    
    # Filter both dataframes to only include common participants
    binary_df = binary_df[binary_df[sample_name].isin(common_participants)]
    values_df = values_df[values_df[sample_name].isin(common_participants)]
    
    # Reset indices to match based on Participant column for easier merging later
    binary_df = binary_df.set_index(sample_name).sort_index()
    values_df = values_df.set_index(sample_name).sort_index()
    
    # Add the Participant column back as a regular column (not index)
    binary_df[sample_name] = binary_df.index
    values_df[sample_name] = values_df.index
    
    # Reset to numeric index while keeping Participant as a column
    binary_df = binary_df.reset_index(drop=True)
    values_df = values_df.reset_index(drop=True)
    
    print(f"After aligning on common participants:")
    print(f"  Binary attributes: {binary_df.shape}")
    print(f"  Values: {values_df.shape}")
    
    # Apply T21 filtering if requested
    if t21_only and t21_column:
        print(f"\nFiltering for T21 patients only...")
        binary_df, values_df = filter_t21_patients(binary_df, values_df, t21_column, sample_name)
    
    return binary_df, values_df

def filter_comorbidities(binary_df, patient_comorbid_threshold, min_comorbids_percent, 
                         max_comorbids_percent, sample_name, remove_comorbidities=None):
    """ This function filters the comorbidity data based on the user defined thresholds. """
    total_samples = binary_df.shape[0]

    if patient_comorbid_threshold < 0:
        raise ValueError("Patient comorbid threshold must be a positive number.")
    if not (0 <= min_comorbids_percent <= 1 and 0 <= max_comorbids_percent <= 1):
        raise ValueError("Comorbid percentage thresholds must be between 0 and 1.")

    # Get all columns except the sample_name column for summing
    comorbidity_cols = [col for col in binary_df.columns if col != sample_name]
    
    # Remove manually specified comorbidities first if provided
    if remove_comorbidities:
        print(f"  Manually removing comorbidities: {remove_comorbidities}")
        cols_before = len(comorbidity_cols)
        # Filter out the specified comorbidities
        comorbidity_cols = [col for col in comorbidity_cols if col not in remove_comorbidities]
        cols_removed = cols_before - len(comorbidity_cols)
        print(f"    Removed {cols_removed} specified comorbidities")
        if cols_removed < len(remove_comorbidities):
            not_found = [c for c in remove_comorbidities if c not in binary_df.columns or c == sample_name]
            if not_found:
                print(f"    Warning: These comorbidities were not found: {not_found}")
    
    # Convert string booleans to actual booleans if needed
    for col in comorbidity_cols:
        if binary_df[col].dtype == 'object':
            binary_df[col] = binary_df[col].map({'True': 1, 'False': 0, True: 1, False: 0})
    
    # Calculate row sums for filtering
    row_sums = binary_df[comorbidity_cols].sum(axis=1)
    filtered_df = binary_df[row_sums >= patient_comorbid_threshold]
    
    print(f"  Patients with >= {patient_comorbid_threshold} comorbidities: {len(filtered_df)} out of {len(binary_df)}")
    
    # Calculate percent occurrence of comorbidities
    df_noname = filtered_df[comorbidity_cols]
    n_comorbid = df_noname.sum(axis=0).to_frame()
    n_comorbid.columns = ["n"]
    n_comorbid["binary_attribute"] = n_comorbid.index
    n_comorbid["percent"] = n_comorbid["n"]/total_samples

    # Filter comorbidities based on prevalence thresholds
    n_comorbid_include = n_comorbid[n_comorbid["percent"]<max_comorbids_percent]
    n_comorbid_include = n_comorbid_include[n_comorbid_include["percent"]>min_comorbids_percent]
    n_comorbid_exclude = n_comorbid[~n_comorbid.index.isin(n_comorbid_include.index)]
    
    # Add manually removed comorbidities to the exclude list for documentation
    if remove_comorbidities:
        for comorbid in remove_comorbidities:
            if comorbid in binary_df.columns and comorbid not in n_comorbid_exclude.index:
                # Add to exclude list with a note that it was manually removed
                manual_row = pd.DataFrame({
                    'n': [0],
                    'binary_attribute': [comorbid],
                    'percent': [-1]  # Use -1 to indicate manual removal
                }, index=[comorbid])
                n_comorbid_exclude = pd.concat([n_comorbid_exclude, manual_row])
    
    # Get the list of comorbidities to keep
    comorbids_to_keep = n_comorbid_include.index.tolist()
    
    # Create filtered dataframe with selected comorbidities AND the Participant column
    columns_to_keep = [sample_name] + comorbids_to_keep
    filtered_binary_df = filtered_df[columns_to_keep]
    
    return filtered_binary_df, n_comorbid_include, n_comorbid_exclude

def filter_gene_expression(values_df, individual_expression_threshold, min_mean_expression, 
                          sample_name, chr21_genes=None, chr21_only=False):
    """ This function filters the gene expression data based on the user defined thresholds. 
    Now includes optional chromosome 21 filtering. """
    
    # Get all columns except the sample_name column
    expression_cols = [col for col in values_df.columns if col != sample_name]
    
    # Apply chromosome filtering FIRST if requested
    if chr21_genes is not None:
        original_count = len(expression_cols)
        expression_cols = filter_genes_by_chromosome(expression_cols, chr21_genes, chr21_only)
        if chr21_only:
            print(f"  Filtered to {len(expression_cols)} chromosome 21 genes (from {original_count} total)")
        else:
            print(f"  Filtered to {len(expression_cols)} non-chromosome 21 genes (from {original_count} total)")
    
    expression_df_nosamples = values_df[expression_cols]

    # Filter individuals (rows) where all gene expression values are below the threshold
    individual_mask = expression_df_nosamples.apply(lambda row: (row >= individual_expression_threshold).any(), axis=1)
    expression_df_filtered = values_df[[sample_name] + expression_cols][individual_mask]
    
    print(f"  Patients with at least one gene >= {individual_expression_threshold}: {len(expression_df_filtered)} out of {len(values_df)}")
    
    # Recalculate after row filtering
    expression_cols = [col for col in expression_df_filtered.columns if col != sample_name]
    expression_df_nosamples = expression_df_filtered[expression_cols]

    # Calculate mean expression
    meansdf = expression_df_nosamples.mean(axis=0).to_frame()
    meansdf.columns = ["mean_value"]
    meansdf["valuename"] = meansdf.index

    if min_mean_expression < 0:
        raise ValueError("Minimum mean expression threshold must be non-negative.")

    # Filter genes based on mean expression threshold
    meansdf_include = meansdf[meansdf["mean_value"] > min_mean_expression]
    meansdf_exclude = meansdf[~meansdf.index.isin(meansdf_include.index)]
    
    # Get the list of genes to keep
    genes_to_keep = meansdf_include.index.tolist()
    
    # Create filtered dataframe with selected genes AND the Participant column
    columns_to_keep = [sample_name] + genes_to_keep
    filtered_values_df = expression_df_filtered[columns_to_keep]

    return filtered_values_df, meansdf_include, meansdf_exclude

def save_filtered_dataframes(filtered_binary_df, 
                            filtered_values_df,
                            n_comorbid_include, 
                            n_comorbid_exclude, 
                            filtered_binary_attribute_file, 
                            filtered_values_file,
                            include_binary_attribute_file, 
                            exclude_binary_attribute_file, 
                            meansdf_include, 
                            meansdf_exclude, 
                            include_values_file, 
                            exclude_values_file):
    """ This function saves the filtered dataframes and include/exclude lists to CSV files. """
    
    # Ensure directories exist
    for filepath in [filtered_binary_attribute_file, filtered_values_file,
                     include_binary_attribute_file, exclude_binary_attribute_file,
                     include_values_file, exclude_values_file]:
        if filepath:  # Check if filepath is not None
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Ensure both dataframes have the same participants before saving
    common_participants = sorted(list(set(filtered_binary_df.index) & set(filtered_values_df.index)))
    if len(common_participants) < len(filtered_binary_df) or len(common_participants) < len(filtered_values_df):
        print(f"Aligning to {len(common_participants)} common participants after filtering")
        filtered_binary_df = filtered_binary_df.loc[common_participants]
        filtered_values_df = filtered_values_df.loc[common_participants]
    
    # Save the filtered dataframes WITH the Participant column
    filtered_binary_df.to_csv(filtered_binary_attribute_file)
    print(f"Filtered binary attribute dataframe saved to {filtered_binary_attribute_file}")
    print(f"  Shape: {filtered_binary_df.shape}")
    print(f"  Columns preview: {list(filtered_binary_df.columns[:5])}..." if len(filtered_binary_df.columns) > 5 else f"  Columns: {list(filtered_binary_df.columns)}")
    
    filtered_values_df.to_csv(filtered_values_file)
    print(f"Filtered values dataframe saved to {filtered_values_file}")
    print(f"  Shape: {filtered_values_df.shape}")
    print(f"  Columns preview: {list(filtered_values_df.columns[:5])}..." if len(filtered_values_df.columns) > 5 else f"  Columns: {list(filtered_values_df.columns)}")
    
    # Verify same number of rows
    assert len(filtered_binary_df) == len(filtered_values_df), "Filtered dataframes must have same number of participants!"
    print(f"  ✓ Both files have {len(filtered_binary_df)} participants")
    
    # Save include/exclude lists
    n_comorbid_include[["binary_attribute"]].to_csv(include_binary_attribute_file, header=False, index=False)
    print(f"Included comorbid list saved to {include_binary_attribute_file}")
    print(f"  Count: {len(n_comorbid_include)}")
    
    n_comorbid_exclude[["binary_attribute"]].to_csv(exclude_binary_attribute_file, header=False, index=False)
    print(f"Excluded comorbid list saved to {exclude_binary_attribute_file}")
    print(f"  Count: {len(n_comorbid_exclude)}")
    
    meansdf_include[["valuename"]].to_csv(include_values_file, header=False, index=False)
    print(f"Included expression list saved to {include_values_file}")
    print(f"  Count: {len(meansdf_include)}")
    
    meansdf_exclude[["valuename"]].to_csv(exclude_values_file, header=False, index=False)
    print(f"Excluded expression list saved to {exclude_values_file}")
    print(f"  Count: {len(meansdf_exclude)}")

def run_filtering(patient_comorbid_threshold, 
                  min_comorbids_percent, 
                  max_comorbids_percent, 
                  individual_expression_threshold, 
                  min_mean_expression, 
                  values_file, 
                  binary_attribute_file, 
                  sample_name, 
                  filtered_values_file,
                  filtered_binary_attribute_file,
                  include_values_file, 
                  include_binary_attribute_file, 
                  exclude_values_file, 
                  exclude_binary_attribute_file,
                  chr21_file=None,
                  chr21_only=False,
                  remove_comorbidities=None,
                  t21_only=False,
                  t21_column='MONDO_complete_trisomy_21'):
    
    """ This is the main function that runs comorbidity and gene expression filtering. """
    try:
        # Load chr21 genes if file provided
        chr21_genes = None
        if chr21_file:
            print(f"\nLoading chromosome 21 genes from {chr21_file}")
            chr21_genes = load_chr21_genes(chr21_file)
        
        # Validate input files and load data - now with T21 filtering
        binary_df, values_df = validate_input_files(values_file, 
                                                    binary_attribute_file, 
                                                    sample_name,
                                                    t21_only,
                                                    t21_column)

        print(f"\nData loaded and aligned:")
        print(f"  Binary attributes shape: {binary_df.shape}")
        print(f"  Values shape: {values_df.shape}")
        if 'Participant' in binary_df.columns:
            print(f"  Participant range: {binary_df['Participant'].min()} to {binary_df['Participant'].max()}")

        # Comorbidity filtering (now with manual removal option)
        print(f"\nFiltering comorbidities...")
        filtered_binary_df, n_comorbid_include, n_comorbid_exclude = filter_comorbidities(
            binary_df, 
            patient_comorbid_threshold, 
            min_comorbids_percent, 
            max_comorbids_percent, 
            sample_name,
            remove_comorbidities)
        
        print(f"Comorbidity filtering results:")
        print(f"  Passing comorbidities: {len(n_comorbid_include)}")
        print(f"  Excluded comorbidities: {len(n_comorbid_exclude)}")
        print(f"  Filtered patient count: {len(filtered_binary_df)}")

        # Gene expression filtering (with chr21 option)
        print(f"\nFiltering gene expression...")
        if chr21_genes is not None:
            filter_type = "chromosome 21 only" if chr21_only else "non-chromosome 21"
            print(f"  Applying {filter_type} filtering")
            
        filtered_values_df, meansdf_include, meansdf_exclude = filter_gene_expression(
            values_df, 
            individual_expression_threshold, 
            min_mean_expression, 
            sample_name,
            chr21_genes,
            chr21_only)
        
        print(f"Gene expression filtering results:")
        print(f"  Passing genes: {len(meansdf_include)}")
        print(f"  Excluded genes: {len(meansdf_exclude)}")
        print(f"  Filtered patient count: {len(filtered_values_df)}")
        
        # Save filtered dataframes and lists
        print(f"\nSaving outputs...")
        save_filtered_dataframes(filtered_binary_df, 
                                 filtered_values_df,
                                 n_comorbid_include, 
                                 n_comorbid_exclude, 
                                 filtered_binary_attribute_file,
                                 filtered_values_file,
                                 include_binary_attribute_file, 
                                 exclude_binary_attribute_file, 
                                 meansdf_include, 
                                 meansdf_exclude, 
                                 include_values_file, 
                                 exclude_values_file)
        
        print("\n✓ Filtering completed successfully!")
    
    except Exception as e:
        print(f"Error occurred during filtering: {e}")
        raise