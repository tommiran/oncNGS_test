import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Script for merging annotated SVs')
parser.add_argument('-i','--input', dest='input', nargs = "+", help='Input annotated SVs')
parser.add_argument('-o','--output', dest='output', help='The output')
args = parser.parse_args()

'''
Merge the annotated SVs into a single table sorted by total supported reads 
'''

# Store tables into list
input_files = []

# Iterate over input files and 
# read as data frames 
for f in args.input:
    tb = pd.read_csv(f, sep = "\t", dtype = str)
    input_files.append(tb)

# Concatenate 
merged_tb = pd.concat(input_files)

# Convert TOTAL_SUPPORTING_READS back to int 
merged_tb["TOTAL_SUPPORTING_READS"] = merged_tb["TOTAL_SUPPORTING_READS"].astype(int) 

# SORT 
merged_tb = merged_tb.sort_values("TOTAL_SUPPORTING_READS", ascending=False)

# Write to output
merged_tb.to_csv(args.output, sep = "\t", index = False)

