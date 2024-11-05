import argparse
import pandas as pd
import numpy as np

opt_parser = argparse.ArgumentParser()
opt_parser.add_argument("-n", "--new_input", dest="new_input", help="Insert a gtf filepath to parse")
opt_parser.add_argument("-s", "--state", dest="state", help="Insert a list with the names of all samples")


options = opt_parser.parse_args()
new_input = options.new_input
state = options.state


#Construct big template dataframe for augmentation
df_temp = pd.read_csv(new_input, sep = "\t", header = 0)
samplenames = list(df_temp.iloc[:,-1])
samplenames = [i for i in samplenames if i != "n"]
name_of_sample = df_temp.columns[-2]
df_construct = pd.DataFrame([[0 for i in range(len(samplenames))] for k in range(df_temp.shape[0])],columns=samplenames)
df_big = pd.concat([df_temp.iloc[:,0:1],df_construct], axis = 1)
df_big.columns = ["Geneid"] + samplenames
for samplename in samplenames:
    if samplename in name_of_sample:
        appended_values = np.array(df_big.loc[:,samplename]) + np.array(df_temp.iloc[:,-2])
        df_big[samplename] = appended_values
        break
df_temp = df_big

try:
    df_temp2 = pd.read_csv(state, sep =  "\t", header = 0)

    df_temp_short = df_temp.iloc[:,1:]
    df_temp2_short = df_temp2.iloc[:,1:]
    df = df_temp_short + df_temp2_short
    df = pd.concat([df_temp["Geneid"],df], axis=1)
    df.to_csv("merged_all_temp.csv", sep = "\t", index = False)
except:
    df_temp.to_csv("merged_all_temp.csv", sep = "\t", index = False)