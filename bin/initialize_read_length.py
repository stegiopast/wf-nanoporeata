import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-s", "--sample_file", dest="sample", help="Insert a sample file to add names to", metavar="FILE")
opt_parser.add_argument("-n", "--sample_id", dest ="sample_id", help="Insert Samplename", metavar="FILE")
opt_parser.add_argument("-m", "--metadata", dest ="metadata", help="Insert metadata file", metavar="FILE")
opt_parser.add_argument("-o", "--output_file", dest="output_file", help="Insert a json file")

options = opt_parser.parse_args()

sample = options.sample
sample_id = options.sample_id
metadata = options.metadata
output_file = options.output_file

line_lengths = []
with open(sample,'r') as file:
    for line in file:
        length = int(line)
        line_lengths.append(length)

metadata_df = pd.read_csv(metadata,header = 0,sep="\t")
samplenames = list(metadata_df["Samples"])

df = pd.DataFrame(np.zeros((40000,len(samplenames))),columns=samplenames)

if sample_id in samplenames:
    print(sample_id)
    for length in line_lengths:
        print(df[sample_id])
        print(length)
        df[sample_id][length] = df[sample_id][length] + 1 

df.to_csv(output_file, sep="\t", index = False)





