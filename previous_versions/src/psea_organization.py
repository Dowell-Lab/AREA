import pandas as pd
import argparse
import time
import itertools


def run_org(outdir,sample_name,values_file, bianary_attribute_file,timestr):
    outfilename=outdir+"organization_"+timestr+".df"
    bianary_attribute_df = pd.read_csv(bianary_attribute_file,index_col=0)
    bianary_attribute_df = bianary_attribute_df.drop(columns=[sample_name])
    values_df = pd.read_csv(values_file, index_col=0)
    values_df = values_df.drop(columns=[sample_name])
    allvalues = values_df.columns.to_list()
    allbiararyattr = bianary_attribute_df.columns.to_list()
    combos_obj = list(itertools.product(allbiararyattr, allvalues))
    tododf = pd.DataFrame(combos_obj)
    tododf.columns = ["binary_attribute", "value"]
    #${outdirname}${bianary_attribute_name}_${value_name}.csv
    tododf["unsorted_rank_file"]=outdir+tododf["binary_attribute"]+"_"+tododf["value"]+".unsoredrank.csv"
    tododf["psea_score_file"]=outdir+tododf["binary_attribute"]+"_"+tododf["value"]+".pseascore.csv"
    print(outfilename)
    tododf.to_csv(outfilename) 
    return outfilename

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run psea on a single set of values and binary attributes')
    parser.add_argument('-baf', '--bianary_attribute_file')
    parser.add_argument('-vf', '--values_file')
    parser.add_argument('-sn', '--sample_name')
    parser.add_argument('-od', '--outdir')
    args = parser.parse_args()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    bianary_attribute_file = args.bianary_attribute_file
    values_file = args.values_file
    outdir = args.outdir
    sample_name = args.sample_name
    run_org(outdir,sample_name,values_file, bianary_attribute_file,timestr)

