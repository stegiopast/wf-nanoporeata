import argparse
import pysam
import pandas as pd
import os

opt_parser = argparse.ArgumentParser()
opt_parser.add_argument("-m", "--metadata", dest="metadata", help="Insert a list with the names of all samples")

options = opt_parser.parse_args()
metadata = options.metadata

metadata_df = pd.read_csv(metadata,header = 0,sep="\t")
samplenames = list(metadata_df["Samples"])

for sample in samplenames:
    bamfile_name = f"{sample}.bam"
    with open(bamfile_name,"w") as file:
        file.write("")

for i in range(1, 97):
    sample = f"{i:02}"
    bamfile_name = f"barcode{sample}.bam"
    with open(bamfile_name,"w") as file:
        file.write("")

with open("unclassified.bam","w") as file:
    file.write("")

