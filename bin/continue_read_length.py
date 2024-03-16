import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-n", "--new_table", dest="new_table", help="Insert a new table", metavar="FILE")
opt_parser.add_argument("-s", "--state", dest="state", help="Insert intermediate state table")

options = opt_parser.parse_args()

new_table = options.new_table
state = options.state

new_table_df = pd.read_csv(new_table,sep="\t",header=0)
state_df = pd.read_csv(state,sep="\t",header=0)

print(new_table)
print(state)

combined_df = new_table_df + state_df
combined_df.to_csv("merged_all_readlengths_temp.csv",sep="\t", index = False)




