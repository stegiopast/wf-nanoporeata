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
print(df_temp)
print(df_temp.columns)
print(df_temp.shape)
df_temp["Name"] = [i.split(sep = "|")[0] for i in df_temp["Name"]]

df = df_temp.iloc[:,0:3]
df.columns = ["Name","Length","EffectiveLength"]
print(df.head())
name_of_sample = input_file
for samplename in samplenames:
    if samplename in name_of_sample:
        print(samplename)
        df[samplename] = np.array(df_temp.iloc[:,-1])
        break
samplenames_transfer_list = samplenames + ["n" for i in range(df.shape[0] - len(samplenames))]
df["all_samplenames"] = samplenames_transfer_list

df.to_csv("merged_all_temp.csv",sep="\t",index = False)
metadata_df.to_csv("metadata.csv",sep="\t",index = False)
