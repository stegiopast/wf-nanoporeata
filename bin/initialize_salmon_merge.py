import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()
opt_parser.add_argument("-i", "--input_file", dest="input_file", help="Insert a gtf filepath to parse")
opt_parser.add_argument("-m", "--metadata", dest="metadata", help="Insert a list with the names of all samples")
opt_parser.add_argument("-o", "--output",dest="output_file", help="Insert a path for the output file", metavar="FILE")

options = opt_parser.parse_args()
input_file = options.input_file
metadata = options.metadata
output_file = options.output_file


metadata_df = pd.read_csv(metadata,header = 0,sep="\t")
samplenames = list(metadata_df["Samples"])
print(options.input_file)
df_temp = pd.read_csv(input_file, sep = "\t", header = 0)
print(df_temp.shape)
df_temp["Name"] = [i.split(sep = "|")[0] for i in df_temp["Name"]]
df_construct = pd.DataFrame([[0 for i in range(len(samplenames))] for k in range(df_temp.shape[0])],columns=samplenames)
print(df_construct.shape)
df = pd.concat([df_temp.iloc[:,0:3],df_construct], axis = 1)
df.columns = ["Name","Length","EffectiveLength"] + samplenames
print(df.head())
name_of_sample = input_file
for samplename in samplenames:
    if samplename in name_of_sample:
        appended_values = np.array(df.loc[:,samplename]) + np.array(df_temp.iloc[:,-1])
        df[samplename] = appended_values
        break
df.to_csv("merged_all_temp.csv",sep="\t",index = False)
metadata_df.to_csv("metadata.csv",sep="\t",index = False)
