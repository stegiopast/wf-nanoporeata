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



if sample_id in samplenames:
    df = pd.DataFrame(np.zeros((30000,1)),columns=[sample_id])
    samplenames_transfer_list = samplenames + ["n" for i in range(df.shape[0] - len(samplenames))]
    df["all_samplenames"] = samplenames_transfer_list
    for length in line_lengths:
        try:
            df[sample_id][length] = df[sample_id][length] + 1 
        except:
            continue
else:
    df = pd.DataFrame(np.zeros((30000,1)),columns=["barcodeXX"])
    samplenames_transfer_list = samplenames + ["n" for i in range(df.shape[0] - len(samplenames))]    
    df["all_samplenames"] = samplenames_transfer_list              

df.to_csv(output_file, sep="\t", index = False)





