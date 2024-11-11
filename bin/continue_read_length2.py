import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-n", "--new_table", dest="new_table", help="Insert a new table", metavar="FILE")
opt_parser.add_argument("-s", "--state", dest="state", help="Insert intermediate state table")

options = opt_parser.parse_args()

new_table = options.new_table
state = options.state

df_temp = pd.read_csv(new_table, sep = "\t", header = 0)

samplenames = list(df_temp.iloc[:,-1])
samplenames = [i for i in samplenames if i != "n"]
print(samplenames)
name_of_sample = df_temp.columns[-2]
print(name_of_sample)

df_big = pd.DataFrame([[0 for i in range(len(samplenames))] for k in range(df_temp.shape[0])],columns=samplenames)
print(df_big)
for samplename in samplenames:
    if samplename in name_of_sample:
        appended_values = np.array(df_big.loc[:,samplename]) + np.array(df_temp.iloc[:,-2])
        df_big[samplename] = appended_values
        break

try:
    new_df = df_big
    state_df = pd.read_csv(state, sep =  "\t", header = 0)    
    print(state_df)
    combined_df = new_df + state_df
    combined_df.to_csv("merged_all_readlengths_temp.csv",sep="\t", index = False)
except:
    df_big.to_csv("merged_all_readlengths_temp.csv",sep="\t", index = False)
    print("New_df")
