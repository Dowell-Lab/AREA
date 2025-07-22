import pandas as pd
import argparse
import numpy as np
from statsmodels.stats.multitest import multipletests


def remove_filenames_from_df(df):
    df = df.drop(columns=["unsorted_rank_file","psea_score_file"])
    return df

def add_adj_Bonferroni(df):
    print (df.shape)
    df['col_num'] = pd.to_numeric(df['pval'], errors='coerce')
    df = remove_filenames_from_df(df)
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
    finaldf = pd.concat([pvalues_numeric_df, pvalues_text_df])
    finaldf = finaldf.drop(columns=["col_num"])
    finaldf = finaldf.sort_values(["p_value_BenjaminiHochberg"])
    return finaldf

def adjust_pvalues(p_values, method):
   return multipletests(p_values, method = method)[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run psea multiple hypothesis correction')
    parser.add_argument('-i', '--input_dataframe_file')
    parser.add_argument('-o', '--output_dataframe_file')
    args = parser.parse_args()
    dataframe_file = args.input_dataframe_file
    dataframe_outfile = args.output_dataframe_file
    df = pd.read_csv(dataframe_file, index_col=0)
    df = add_adj_Bonferroni(df)
    df.to_csv(dataframe_outfile)

