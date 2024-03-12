#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

file1=args[1]
file2=args[2]
DT1=read.table(file1,sep="\t",header=T)
DT2=read.table(file2,sep="\t",header=T)

print(dim(DT1))
print(dim(DT2))

output <- as.numeric(DT1[,7]) + as.numeric(DT2[,7])
print(output)
sum_df <- data.frame(DT1[,1:6], Output = output)


write.table(sum_df,"feature_counts_latest.csv", row.names = FALSE,sep="\t")
