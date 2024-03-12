import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()
opt_parser.add_argument("-n", "--new_input", dest="new_input", help="Insert a gtf filepath to parse")
opt_parser.add_argument("-s", "--state", dest="state", help="Insert a list with the names of all samples")
#opt_parser.add_argument("-m", "--metadata", dest="metadata", help="Insert a metadata sheet")


options = opt_parser.parse_args()
new_input = options.new_input
state = options.state


df_temp = pd.read_csv(new_input, sep = "\t", header = 0)
df_temp2 = pd.read_csv(state, sep =  "\t", header = 0)

df_temp_short = df_temp.iloc[:,3:]
df_temp2_short = df_temp2.iloc[:,3:]
print(df_temp.head())
print(df_temp2.head())

df = df_temp_short + df_temp2_short
df = pd.concat([df_temp[["Name","Length","EffectiveLength"]],df], axis=1)
df.to_csv("merged_all_temp.csv", sep = "\t", index = False)