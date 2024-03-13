import argparse
import pandas as pd
import numpy as np
import os



opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-s", "--sample_file", dest="sample", help="Insert a sample file to add names to", metavar="FILE")
opt_parser.add_argument("-d", "--percentages_output_path",dest="percentages", help="Insert an output path for percenatges", metavar="FILE")
opt_parser.add_argument("-o", "--output_path",dest="output", help="Insert a template file to extract names from", metavar="FILE")


options = opt_parser.parse_args()

sample = options.sample
output_path = options.output
p_output_path = options.percentages

#Read FeatureCount count table 
sample_df = pd.read_csv(sample,header = 0, index_col = 0, sep = "\t")

#Check for necessary paths
if os.path.getsize(output_path) < 3:
    output_df = pd.DataFrame()
else:
    output_df = pd.read_csv(output_path,header=0, sep = "\t")

if os.path.getsize(p_output_path) < 3:
    p_output_df = pd.DataFrame()
else:
    p_output_df = pd.read_csv(p_output_path, header=0, sep = "\t")


#Read the relative read counts for every sample in the FeatureCount table
samplenames = []
columns_of_percentages = [] 
for i in range(len(sample_df.iloc[0,:])):
    column = sample_df.iloc[:,i]
    if not type(column[0]) == type(""):
        name = column.name
        samplenames.append(name)
        sum = column.sum()
        percentages = []
        for i in column:
            if sum == 0:
                percentages.append(0)
            else:
                percentages.append(float(i)/sum)
        columns_of_percentages.append(percentages)

#Append relative read counts to a new dataframe
new_percentages_row = pd.DataFrame(columns = samplenames, index = [0])
for i,val in enumerate(samplenames):
    new_percentages_row.loc[0,val] = columns_of_percentages[i]
p_output_df = pd.concat([p_output_df,new_percentages_row])
p_output_df = p_output_df.reset_index(drop = True)
p_output_df.to_csv(p_output_path, sep="\t", index = 0)

#Compare relative read counts of this iteration to the last iteration
if len(p_output_df) > 1:
    mean_list = []
    for i,val in enumerate(samplenames):
        last = list(p_output_df[val])[-1]
        before_last = list(p_output_df[val])[-2].replace("[","")
        before_last = before_last.replace("]","")
        before_last = before_last.replace(",","")
        before_last = before_last.split()
        b_last = []
        for i in before_last:
            b_last.append(float(i))
        tmp_list = []
        for k in range(len(last)):
            difference = abs(abs(last[k]) - abs(b_last[k]))
            tmp_list.append(difference)
        tmp_mean = np.mean(tmp_list)
        mean_list.append(tmp_mean)
    print(mean_list)

try:
    new_mean_df = pd.DataFrame([mean_list], columns = samplenames)
except:
    new_mean_df = pd.DataFrame([[0 for i in samplenames]], columns = samplenames)

#Write output
output_df = pd.concat([output_df,new_mean_df])
print(output_df)
output_df.to_csv(output_path, sep="\t", index = 0)



