import argparse
import pandas as pd
import gtfparse as gtf_parse

opt_parser = argparse.ArgumentParser()
opt_parser.add_argument("-i", "--input", dest="input_file", help="Insert a gtf filepath to parse", metavar="FILE")
opt_parser.add_argument("-o", "--output",dest="output_file", help="Insert a path for the output file", metavar="FILE")
options = opt_parser.parse_args()
input_file = options.input_file
output_file = options.output_file

"""Convert GTF file into csv format"""
print(options.input_file)
df = gtf_parse.read_gtf(input_file)
df = df.to_pandas()
df = df[['gene_id', 'gene_name', 'transcript_id', 'transcript_name']]
df.to_csv(output_file)
