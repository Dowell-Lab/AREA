import argparse
import numpy as np
import pandas as pd
import psea_organization
import make_rankable_file
import psea_apval
import psea_core
from multiprocessing import Pool
import time

#def add(row):
#   return row[0]+row[1]+row[2]

#df['new_col'] = df.apply(add, axis=1)


def several_rank_and_cores(dfs_and_args):
    df,sample_name, bianary_attribute_file, values_file = dfs_and_args
    collection_onegeneNES = []
    collection_onegenep = []
    LOG_EVERY_N=1000
    for index, row in df.iterrows():
        if (index % LOG_EVERY_N) == 0:
             print (index)
        shouldIrun=row["runpsea"]
        if shouldIrun=="included":
            outputfile=row["psea_score_file"]
            value_name=row["value"]
            bianary_attribute_name=row["binary_attribute"]	
            gene_cormorbid_data = make_rankable_file.run_merge(outputfile,sample_name,values_file, bianary_attribute_file, value_name, bianary_attribute_name, low_memory=False)
            #now I need to adapt make_rankable_file.run_merge run off a dataframe
            onegeneNES, onegenep = psea_core.psea_heart(gene_cormorbid_data, outputfile)
            collection_onegeneNES.append(onegeneNES)
            collection_onegenep.append(onegenep)
        else:
            collection_onegeneNES.append("NA")
            collection_onegenep.append("NA")
    df["NES"]=collection_onegeneNES
    df["pval"]=collection_onegenep
    print("final dfs", df)
    return df
        

def change_master_runlist_via_include_exclude(df, include_values_file, include_bianary_attribute_file, exclude_values_file, exclude_bianary_attribute_file):
     print ("before include and exclude files", df["runpsea"].value_counts())
     include_values_list = []
     exclude_values_list = []
     include_ba_list = []
     exclude_ba_list = []
     if include_values_file!=False:
        include_values_df = pd.read_csv(include_values_file, names=["value"])
        include_values_list = include_values_list+include_values_df["value"].to_list()
     if include_bianary_attribute_file!=False:
        include_ba_df = pd.read_csv(include_bianary_attribute_file, names=["binary_attribute"])
        include_ba_list = include_ba_list+include_ba_df["binary_attribute"].to_list()
     if exclude_values_file!=False:
        exclude_values_df = pd.read_csv(exluded_values_file, names=["value"])
        exclude_values_list = exclude_values_list+exclude_values_df["value"].to_list()
     if exclude_bianary_attribute_file!=False:
        exlude_ba_df = pd.read_csv(exclude_bianary_attribute_file, names=["binary_attribute"])
        exclude_ba_list = exclude_ba_list+exlude_ba_df["binary_attribute"].to_list()
     include_lens = len(include_values_list)+len(include_ba_list)
     exclude_lens = len(exclude_values_list)+len(exclude_ba_list)
     if include_lens>0 and exclude_lens>0:
         print ("You can only use include or exclude lists. Nothing will be run becuase both include and exclude list are used.")
         df["runpsea"]=="exclude_all" 
     elif include_lens>0:
         if len(include_values_list)>0 and len(include_ba_list)>0:
             df["runpsea"] = np.where((df["value"].isin(include_values_list) & df["binary_attribute"].isin(include_ba_list)), 'included','excluded')
             print ("here1")
         elif len(include_values_list)>0:
             df["runpsea"] = np.where(df["value"].isin(include_values_list), 'included','excluded')
             print ("here2")
         elif len(include_ba_list)>0:
             df["runpsea"] = np.where(df["binary_attribute"].isin(include_ba_list), 'included','excluded')
             print ("here3")
     elif exclude_lens>0:
         df["runpsea"] = np.where((df["value"].isin(exclude_values_list) | df["binary_attribute"].isin(exculde_ba_list)), 'included','excluded')
         print ("here4")
     print ("after include and exclude files", df["runpsea"].value_counts())
     return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run psea on a set of values and binary attributes')
    parser.add_argument('-baf', '--bianary_attribute_file')
    parser.add_argument('-vf', '--values_file')
    parser.add_argument('-sn', '--sample_name')
    parser.add_argument('-od', '--outdir')
    parser.add_argument('-p', '--processes', default=4,type=int)
    parser.add_argument('-ivf', '--include_values_file', default=False)
    parser.add_argument('-ibaf', '--include_bianary_attribute_file', default=False)
    parser.add_argument('-evf', '--exclude_values_file', default=False)
    parser.add_argument('-ebaf', '--exclude_bianary_attribute_file', default=False)
    args = parser.parse_args()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    bianary_attribute_file = args.bianary_attribute_file
    values_file = args.values_file
    outdir = args.outdir
    sample_name = args.sample_name
    master_org_filename = psea_organization.run_org(outdir,sample_name,values_file, bianary_attribute_file,timestr)
    master_org_df = pd.read_csv(master_org_filename, index_col=0)
    master_org_df["runpsea"] = "included"
    print (master_org_df)
    if args.include_values_file!=False or args.include_bianary_attribute_file!=False or args.exclude_values_file!=False or args.exclude_bianary_attribute_file!=False:
        master_org_df = change_master_runlist_via_include_exclude(master_org_df,args.include_values_file, args.include_bianary_attribute_file, args.exclude_values_file, args.exclude_bianary_attribute_file) 
    n_processes = args.processes
    df_list = np.array_split(master_org_df, n_processes) # split the dataframe into 161 separate tables
    dfs_and_args = [(df, sample_name, bianary_attribute_file, values_file) for df in df_list] 
    with Pool(processes=n_processes) as p:
        new_df_list = p.map(several_rank_and_cores, dfs_and_args)    
    outfilename=outdir+"psea_scores_"+timestr
    print(outfilename)
    finalresult = pd.concat(new_df_list)
    finalresult.to_csv(outfilename+".csv")
    dataframe_outfile =  outfilename+".adjpval.csv"
    finalresult = psea_apval.add_adj_Bonferroni(finalresult)
    finalresult.to_csv(dataframe_outfile) 
    wf = open(outfilename+".info", "w")
    line = "this will be where I put all the arguments used when running this"
    wf.write(line)
    wf.close()
