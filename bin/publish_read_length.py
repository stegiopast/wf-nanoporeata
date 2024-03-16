import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-i", "--input_table", dest="input_table", help="Insert a new table", metavar="FILE")

options = opt_parser.parse_args()

input_table = options.input_table

df = pd.read_csv(input_table,sep="\t",header=0)

samplenames = df.columns

for sample in samplenames:
    if df[sample].sum() > 0:
        with open(f"{sample}_read_length_pass.txt","w") as file:
            # file.write("Length\n")
            # for index,value in enumerate(df[sample]):
            #     for i in range(int(value)):
            #         file.write(f"{index}\n")
            [[file.write(f"{index}\n") for i in range(int(value))] for index,value in enumerate(df[sample])]
