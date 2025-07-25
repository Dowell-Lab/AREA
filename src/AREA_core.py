import argparse
import pandas as pd
import importlib
import sys
import itertools
import os
from multiprocessing import Pool
from statsmodels.stats.multitest import multipletests
import scipy
import time 
from scipy import stats

# Try to import prefilter functions if available
try:
    from filter_functions import run_filtering
    PREFILTER_AVAILABLE = True
except ImportError:
    PREFILTER_AVAILABLE = False
    print("Prefilter module not found. Prefiltering features disabled.")

#try: 
#    import cupy as cp
#    gpu_available = True
#except ImportError:
#    import numpy as np
#    gpu_available = False
#    print ("gpu not available, using CPU")

# Set np or cp dynamically based on GPU availability.
#array_lib = cp if gpu_available else np

try:
    import cupy as cp
    gpu_available = True
    array_lib = cp
    print("GPU available, using cupy")
except ImportError:
    import numpy as np
    gpu_available = False
    array_lib = np
    print("GPU not available, using CPU")

def area_score_norm(ar_ticks, verbose=False):
    '''Calculate the GSEA like enrichment score using the comorbidity
    occurance in the rank as our set'''

    ar_score = array_lib.array([1 if i > 0 else 0 for i in ar_ticks])  # Use array_lib (either cp or np)
    total = float(array_lib.sum(ar_score))  # Use array_lib sum

    binwidth = 1.0 / float(len(ar_score))
    normalized_ar_score = array_lib.multiply(array_lib.divide(ar_score, total), binwidth)  # Normalization with array_lib

    cumscore = array_lib.cumsum(normalized_ar_score)  # Cumulative sum with array_lib
    trend = array_lib.append(array_lib.arange(0, 1, 1.0 / float(len(cumscore) - 1)), 1.0)
    trend = array_lib.multiply(trend, binwidth)  # Multiplying trend by binwidth

    actual_es = (array_lib.trapz(cumscore) - array_lib.trapz(trend)) * 2

    if verbose:
        print("Binwidth                        :"+str(binwidth))
        print("Sum of normalized binary scores :"+str(array_lib.sum(normalized_ar_score)))
        print("Sum cumulative score            :"+str(array_lib.sum(cumscore)))
        print("Sum of trend                    :"+str(array_lib.sum(trend)))
        print("Len of trend                    :"+str(len(trend)))
        print("Actual Area                     :"+str(actual_es))

    return [actual_es, normalized_ar_score, trend, cumscore]

def permute_area_norm(ar_ticks, permutations=1000, seed=42):
    '''Generates permutations of the ranks and calculates AUC for each
        permutation.'''

    array_lib.random.seed(seed=seed)
    ar_ticks = array_lib.array(ar_ticks)
    es_permute = []
    for i in range(permutations):
        random_ar_score = array_lib.random.permutation(ar_ticks)  # Array_lib permutation
        es_permute_score = area_score_norm(random_ar_score)  # Use area_score_norm
        es_permute.append(es_permute_score[0])

    return es_permute

def calculateNESpval(actualES, simES):
    '''Calculate NES and p-value'''
    if actualES > 0:
            simESsubset = array_lib.array([x for x in simES if x < 0])
            mu = array_lib.mean(simESsubset)  # Use array_lib mean
            NES = actualES / mu
            sigma = array_lib.std(simESsubset)  # Use array_lib std
            p = 1 - scipy.stats.norm.cdf(actualES, mu, sigma)
    else:
            simESsubset = array_lib.array([x for x in simES if x > 0])
            mu = array_lib.mean(simESsubset)  # Use array_lib mean
            NES = -(actualES / mu)
            sigma = array_lib.std(simESsubset)  # Use array_lib std
            p = scipy.stats.norm.cdf(actualES, mu, sigma)
    return NES, p

def is_unique(s):
    a = s.to_numpy() # s.values (pandas<0.24)
    return (a[0] == a).all()

def pre_organize_run(outdir, orgfile, commoncolumn, value_file, binary_attribute_file, keepbas_file=None, keepvalues_file=None, ignorebas_file=None, ignorevalues_file=None,verbose=False):
    #read the files
    badf = pd.read_csv(binary_attribute_file, index_col=0)
    rankorderdf = pd.read_csv(value_file, index_col=0)

    #gather all pairs
    valuelabels = [colname for colname in rankorderdf.columns if colname!=commoncolumn]
    balabels = [colname for colname in badf.columns if colname!=commoncolumn]
    product = list(itertools.product(valuelabels, balabels))
    # Create DataFrame from the product
    torundf = pd.DataFrame(product, columns=['valuelabel', 'balabel'])
    torundf["plan"] = "run_area"
    if isinstance(ignorebas_file, str):
        if verbose:
            print ("excluding attributes from", ignorebas_file)
        ignoreba = pd.read_csv(ignorebas_file, names=["ignoreba"])
        torundf.loc[torundf['balabel'].isin(ignoreba["ignoreba"]), 'plan'] = 'ignored due to user input ignore file'
    if isinstance(ignorevalues_file, str):    
        if verbose:
            print ("excluding values from", ignorevalues_file)
        ignorevalue = pd.read_csv(ignorevalues_file, names=["ignorevalue"])
        torundf.loc[torundf['valuelabel'].isin(ignorevalue["ignorevalue"]), 'plan'] = 'ignored due to user input ignore file'
    if isinstance(keepbas_file, str):    
        keepba = pd.read_csv(keepbas_file, names=["keepba"])
        if verbose:
            print ("excluding attributes not in ", keepbas_file)
            print(keepba["keepba"])
        torundf.loc[~torundf['balabel'].isin(keepba["keepba"].to_list()), 'plan'] = 'ignored due to user input not in keep file'
    if isinstance(keepvalues_file, str):   
        keepvalue = pd.read_csv(keepvalues_file, names=["keepvalue"])
        if verbose:
            print ("excluding values not in ", keepvalues_file)
            print(keepvalue["keepvalue"])
        torundf.loc[~torundf['valuelabel'].isin(keepvalue["keepvalue"].to_list()), 'plan'] = 'ignored due to user input not in keep file'
    for colname in valuelabels:
        if is_unique(rankorderdf[colname])==True:
            torundf.loc[torundf['valuelabel']==colname, 'plan'] = 'excluded_by_area all values identical'
    for colname in balabels:
        if is_unique(badf[colname])==True:
            torundf.loc[torundf['balabel']==colname, 'plan'] = 'excluded_by_area all binary attributes identical'
    torundf.to_csv(outdir+orgfile)


def process_balabel(balabel, todolist, outdir, commoncolumn, value_file, binary_attribute_file):
    '''Process a single balabel and return the results'''
    todolistthisba = todolist[todolist["balabel"] == balabel]
    valuecols = todolistthisba["valuelabel"].to_list()
    baNESpvals = run_a_ba(balabel, valuecols, outdir, commoncolumn, value_file, binary_attribute_file)
    return baNESpvals

def org_to_pval(outdir, orgfile, commoncolumn, value_file, binary_attribute_file, n_processes=4):
    todolist = pd.read_csv(outdir + orgfile)
    # Remove things we are not going to run because of ignores or being non-unique
    todolist = todolist[todolist["plan"] == "run_area"]
    print(todolist["balabel"].unique())
    # Set up multiprocessing
    with Pool(processes=n_processes) as pool:
        # Use map to apply the `process_balabel` function in parallel
        results = pool.starmap(process_balabel, [(balabel, todolist, outdir, commoncolumn, value_file, binary_attribute_file) for balabel in todolist["balabel"].unique()])
    
    # Collect the results into one DataFrame
    collectdfs = pd.concat(results, ignore_index=True)
    
    # Optionally, save or return the collected DataFrame
    collectdfs.to_csv(outdir + orgfile+'beforeadjpval.csv', index=False)


def pval_to_adjpvals(outdir, orgfile):
    df = pd.read_csv(outdir + orgfile+'beforeadjpval.csv')
    finaldf = add_adj_pvals(df)
    finaldf.to_csv(outdir + orgfile+'.adjpval.csv', index=False)


def run_a_ba(ba,valuecols, outdir, commoncolumn, value_file, binary_attribute_file):
    badf = pd.read_csv(binary_attribute_file, index_col=0)
    rankorderdf = pd.read_csv(value_file, index_col=0)
    bacols_common = [ba, commoncolumn]
    valuecols_common = valuecols+[commoncolumn]
    subbadf = badf[bacols_common]
    subrankorderdf = rankorderdf[valuecols_common]
    mergedf = subbadf.merge(subrankorderdf, on=commoncolumn)
    df = mergedf.sample(frac=1).reset_index(drop=True)
    thiscomorbidity_binary = df[ba].to_list()
    simES_norm = permute_area_norm(thiscomorbidity_binary, permutations=1000)
    lines = []
    for colname in valuecols:
        df = df.sort_values(colname)
        thiscomorbidity_binary = df[ba].to_list()
        actualES_norm, normalized_pc_score_norm, thistrend_norm, thiscumscore_norm = area_score_norm(thiscomorbidity_binary)
        onegeneNES, onegenepval = calculateNESpval(actualES_norm, simES_norm)
        line = [colname, onegeneNES, onegenepval]
        lines.append(line)
    baNESpvals = pd.DataFrame(lines, columns = ["value","NES", "pval"])
    baNESpvals["binary_attribute"]=ba
    baNESpvals = baNESpvals[["binary_attribute", "value","NES", "pval"]]
    return baNESpvals


def add_adj_pvals(df):
    df['col_num'] = pd.to_numeric(df['pval'], errors='coerce')
    #split to two data frames one with numeric pvalues and the other with text
    pvalues_numeric_df=df[~df["col_num"].isna()].copy()
    pvalues_text_df=df[df["col_num"].isna()].copy()
    pvalues_text_df['pval']=pd.NA
    pvalues_numeric_df['pval'] = pvalues_numeric_df["col_num"]
    n_of_psea = pvalues_numeric_df.shape[0]
    #add adjusted pvalues to numeric pvalues and set non-numeric to NA
    pvalues_numeric_df['p_value_bonf'] = adjust_pvalues(pvalues_numeric_df['pval'], 'bonferroni')
    pvalues_text_df['p_value_bonf'] =pd.NA
    pvalues_numeric_df['p_value_holm'] = adjust_pvalues(pvalues_numeric_df['pval'], 'holm')
    pvalues_text_df['p_value_holm']= pd.NA
    pvalues_numeric_df['p_value_BenjaminiHochberg'] = adjust_pvalues(pvalues_numeric_df['pval'], 'fdr_bh')
    pvalues_text_df['p_value_BenjaminiHochberg']= pd.NA
    pvalues_numeric_df['p_value_BenjaminiYekutieli'] = adjust_pvalues(pvalues_numeric_df['pval'], 'fdr_by')
    pvalues_text_df['p_value_BenjaminiYekutieli']= pd.NA
    #concat the files back to gether
    finaldf = pvalues_numeric_df
    finaldf = finaldf.drop(columns=["col_num"])
    finaldf = finaldf.sort_values(["p_value_BenjaminiHochberg"])
    return finaldf

def adjust_pvalues(p_values, method):
   return multipletests(p_values, method = method)[1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run AREA analysis with optional prefiltering')
    
    # Core AREA arguments
    parser.add_argument('-baf', '--bianary_attribute_file', required=True)
    parser.add_argument('-vf', '--values_file', required=True)
    parser.add_argument('-cc', '--common_column_name', required=True)
    parser.add_argument('-od', '--outdir', required=True)
    parser.add_argument('-p', '--processes', default=4, type=int)
    parser.add_argument('-ivf', '--include_values_file', default=False)
    parser.add_argument('-ibaf', '--include_bianary_attribute_file', default=False)
    parser.add_argument('-evf', '--exclude_values_file', default=False)
    parser.add_argument('-ebaf', '--exclude_bianary_attribute_file', default=False)
    
    # Optional prefiltering arguments (only available if prefilter module exists)
    if PREFILTER_AVAILABLE:
        parser.add_argument('--enable_prefilter', action='store_true', 
                            help='Enable automatic prefiltering of data')
        parser.add_argument('--patient_comorbid_threshold', type=int, default=1, 
                            help='Minimum number of comorbidities a patient should have (default: 1)')
        parser.add_argument('--min_comorbids_percent', type=float, default=0.05, 
                            help='Minimum percent prevalence of comorbidity to retain (default: 0.05)')
        parser.add_argument('--max_comorbids_percent', type=float, default=0.95, 
                            help='Maximum percent prevalence of comorbidity to retain (default: 0.95)')
        parser.add_argument('--min_mean_expression', type=float, default=1.0, 
                            help='Minimum mean gene expression to retain gene (default: 1.0)')
        parser.add_argument('--individual_expression_threshold', type=int, default=10, 
                            help='Minimum individual expression threshold (default: 10)')
    
    parser.add_argument('--verbose', action='store_true', 
                        help='Enable verbose output')
    
    args = parser.parse_args()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    orgfile = "area_scores_" + timestr
    
    # Extract arguments
    binary_attribute_file = args.bianary_attribute_file
    value_file = args.values_file
    outdir = args.outdir
    commoncolumn = args.common_column_name
    n_processes = args.processes
    
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # Handle include/exclude files
    keepvalues_file = args.include_values_file if args.include_values_file != False else None
    keepbas_file = args.include_bianary_attribute_file if args.include_bianary_attribute_file != False else None
    ignorebas_file = args.exclude_bianary_attribute_file if args.exclude_bianary_attribute_file != False else None
    ignorevalues_file = args.exclude_values_file if args.exclude_values_file != False else None
    
    # Run prefiltering if enabled and available
    if PREFILTER_AVAILABLE and hasattr(args, 'enable_prefilter') and args.enable_prefilter:
        print("=== PREFILTERING STAGE ===")
        
        # Generate filter file names
        prefilter_dir = os.path.join(outdir, "prefilter_" + timestr + "/")
        include_values_file = os.path.join(prefilter_dir, "include_values.csv")
        include_binary_attribute_file = os.path.join(prefilter_dir, "include_binary_attributes.csv")
        exclude_values_file = os.path.join(prefilter_dir, "exclude_values.csv")
        exclude_binary_attribute_file = os.path.join(prefilter_dir, "exclude_binary_attributes.csv")
        
        try:
            # Run prefiltering using imported function
            run_filtering(
                patient_comorbid_threshold=args.patient_comorbid_threshold,
                min_comorbids_percent=args.min_comorbids_percent,
                max_comorbids_percent=args.max_comorbids_percent,
                individual_expression_threshold=args.individual_expression_threshold,
                min_mean_expression=args.min_mean_expression,
                values_file=value_file,
                binary_attribute_file=binary_attribute_file,
                sample_name=commoncolumn,
                include_values_file=include_values_file,
                include_binary_attribute_file=include_binary_attribute_file,
                exclude_values_file=exclude_values_file,
                exclude_binary_attribute_file=exclude_binary_attribute_file
            )
            
            # Use prefiltered files for subsequent analysis
            keepvalues_file = include_values_file
            keepbas_file = include_binary_attribute_file
            print(f"Prefiltering completed. Using filtered files:")
            print(f"  Values: {keepvalues_file}")
            print(f"  Binary attributes: {keepbas_file}")
            
        except Exception as e:
            print(f"Prefiltering failed: {e}")
            print("Continuing without prefiltering.")
    
    elif hasattr(args, 'enable_prefilter') and args.enable_prefilter and not PREFILTER_AVAILABLE:
        print("Warning: Prefiltering requested but filter_functions module not available.")
        print("Please ensure filter_functions.py is in the same directory.")
        print("Continuing without prefiltering.")
    
    # Run main AREA analysis
    print("=== AREA ANALYSIS STAGE ===")
    
    # Fix the argument passing bug from original
    pre_organize_run(outdir, orgfile, commoncolumn, value_file, binary_attribute_file, 
                     keepbas_file=keepbas_file, keepvalues_file=keepvalues_file, 
                     ignorebas_file=ignorebas_file, ignorevalues_file=ignorevalues_file, 
                     verbose=args.verbose)
    
    org_to_pval(outdir, orgfile, commoncolumn, value_file, binary_attribute_file, n_processes=n_processes)
    pval_to_adjpvals(outdir, orgfile)
    
    print(f"AREA analysis completed. Results saved with prefix: {orgfile}")