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
construct_df = pd.DataFrame([[samplenames[i],0,0] for i in range(len(samplenames))], columns = ["Sample","num_reads","num_mapped_reads"])
print(construct_df)
input_df = pd.read_csv(input_file, sep=";", header = 0)
print(input_df)
construct_df.loc[construct_df["Sample"] == input_df["Sample"][0],"num_reads"] += input_df["num_reads"][0]
construct_df.loc[construct_df["Sample"] == input_df["Sample"][0],"num_mapped_reads"] += input_df["num_mapped_reads"][0]
construct_df.to_csv(output_file,sep="\t",index=False)
