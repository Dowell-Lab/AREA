import pandas as pd
import argparse

def run_merge(outputfile,sample_name,values_file, bianary_attribute_file, value_name, bianary_attribute_name, low_memory=True):
    bianary_attribute_df = pd.read_csv(bianary_attribute_file)
    values_df = pd.read_csv(values_file)
    bianary_attribute_df = bianary_attribute_df[[sample_name, bianary_attribute_name]]
    values_df = values_df[[sample_name, value_name]]
    df = bianary_attribute_df.merge(values_df, how='inner', on=sample_name)
    df=df.rename(columns={sample_name:"sample", value_name:"value", bianary_attribute_name:'binary_attribute'})
#,rank_gene,sample,value,binary_attribute
    if low_memory==True:
        df.to_csv(outputfile)
    else:
        return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run psea on a single set of values and binary attributes')
    parser.add_argument('-baf', '--bianary_attribute_file')
    parser.add_argument('-ban', '--bianary_attribute_name')
    parser.add_argument('-vf', '--values_file')
    parser.add_argument('-vn', '--value_name')
    parser.add_argument('-of', '--outputfile')
    parser.add_argument('-sn', '--sample_name')
    args = parser.parse_args()
    bianary_attribute_file = args.bianary_attribute_file
    bianary_attribute_name = args.bianary_attribute_name
    values_file = args.values_file
    value_name = args.value_name
    outputfile = args.outputfile
    sample_name = args.sample_name
    print (bianary_attribute_file, bianary_attribute_name, values_file, value_name, outputfile, sample_name)
    run_merge(outputfile,sample_name,values_file, bianary_attribute_file, value_name, bianary_attribute_name)


