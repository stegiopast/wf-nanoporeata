#!/usr/bin/env bash


###### Function to extract read length from fastq files 
fastq_file="$1" # $input_dir/*/*/*/*.fastq.gz
outFile="$2" # ERRxyz

#Filter read length information from fastq files
if [[ "$fastq_file" == *".gz" ]]
then 
    zcat $fastq_file | awk 'NR%4==2' | awk '{ print length }' >> $outFile
else
    cat $fastq_file | awk 'NR%4==2' | awk '{ print length }' >> $outFile
fi
