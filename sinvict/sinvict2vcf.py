import argparse
import pandas as pd
import os 
from pyfasta import Fasta

parser = argparse.ArgumentParser(description='Script for converting sinvict output files to a VCF file')
parser.add_argument('-i','--input', dest='input', nargs = "+", help='Input sinvict file')
parser.add_argument('-r','--reference', dest='reference', help='Reference fasta')
parser.add_argument('-o','--output', dest='output', help='Output VCF file')
args = parser.parse_args()

'''
Columns : 
0. Chromosome Name
1. Position
2. Sample Name (readcount file name)
3. Reference Base
4. Read Depth
5. Mutated Base(s)
6. Number of reads supporting the mutation
7. Variant Allele Percentage
8. Number of reads mapped to the + strand (in format "+:value")
9. Number of reads mapped to the - strand (in format "-: value")
10. Average position of the base on reads as a fraction (0.00 - 1.00)
11. "Somatic" or "Germline", predicted based on variant allele frequency (do not take as the ground truth)
'''

'''
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
'''

'''
Match chromosome with fasta header 
'''
def matchChromToHeader(chrom, header_keys):
    for key in header_keys:
        if chrom == key.split()[0]:
            return(key)
        else:
            pass 


# Read in fasta 
fasta = Fasta(args.reference)

# Read in sinvict files
header = ["CHROM", "POS", "ID", "REF", "ALT", "INFO", "FORMAT", "SAMPLE"]

# Prepare variant dictionaries to determine the filtering status 
dc_variants = {"low_dp" : [], 
              "strand_bias" : [],
              "pos_bias" : [],
              "signal_to_noise_ratio" : [],
              "homopolymer" : []}

# Store the unfiltered variants into a list 
all_variants = []

for f in args.input:
    if "level1" in f:
        with open(f) as infile:
            for ln in infile:
                all_variants.append(ln)

    elif "level2" in f:
        if os.stat(f).st_size > 0:
            with open(f) as infile:
                for ln in infile:
                    cols = ln[:-1].split("\t")
                    varid = ":".join([cols[0],cols[1],cols[3], cols[5]])
                    dc_variants["low_dp"].append(varid)

    elif "level3" in f:
        if os.stat(f).st_size > 0:
            with open(f) as infile:
                for ln in infile:
                    cols = ln[:-1].split("\t")
                    varid = ":".join([cols[0],cols[1],cols[3], cols[5]])
                    dc_variants["strand_bias"].append(varid)

    elif "level4" in f:
        if os.stat(f).st_size > 0:
            with open(f) as infile:
                for ln in infile:
                    cols = ln[:-1].split("\t")
                    varid = ":".join([cols[0],cols[1],cols[3], cols[5]])
                    dc_variants["pos_bias"].append(varid)

    elif "level5" in f:
        if os.stat(f).st_size > 0:
            with open(f) as infile:
                for ln in infile:
                    cols = ln[:-1].split("\t")
                    varid = ":".join([cols[0],cols[1],cols[3], cols[5]])
                    dc_variants["signal_to_noise_ratio"].append(varid)

    elif "level6" in f:
        if os.stat(f).st_size > 0:
            with open(f) as infile:
                for ln in infile:
                    cols = ln[:-1].split("\t")
                    varid = ":".join([cols[0],cols[1],cols[3], cols[5]])
                    dc_variants["homopolymer"].append(varid)


# Normalise the variants according to VCF standards 
def normaliseVariant(chrom, pos, ref, alt, fasta):
    # Find sign 
    if "-" in alt:
        # Deletion
        pos = int(pos) - 1
        # Find the reference base at pos 
        header_match = matchChromToHeader(chrom, list(fasta.keys()))
        # Alt 
        alt_final = fasta[header_match][pos - 1]
        # Find the reference  
        ref_final = alt_final + alt[1:] 
        return(str(pos), ref_final, alt_final)
    elif "+" in alt:
        return(pos, ref, ref + alt[1:])
    else:
        return(pos, ref, alt)

# FIXED format column
format_col = "GT:DP:AD:AF"

# Output 
out = open(args.output, "w")

# Iterate over all variants and check which filters are passed 
for ln in all_variants:

    # Split columns 
    cols = ln[:-1].split("\t")
    
    # Extract chromsome, position, reference and alternative alleles
    chrom,pos,ref,alt = cols[0],cols[1],cols[3],cols[5]

    # Construct the varid 
    varid = ":".join([chrom, pos, ref, alt])

    # Normalise variant 
    pos, ref, alt = normaliseVariant(chrom, pos, ref, alt, fasta)

    # Determine the filtering status
    if varid in dc_variants["homopolymer"]:
        filter_col = "PASS"
    elif varid in dc_variants["signal_to_noise_ratio"]:
        filter_col = "Homopolymer"
    elif varid in dc_variants["pos_bias"]:
        filter_col = "signal_to_noise_ratio"
    elif varid in dc_variants["strand_bias"]:
        filter_col = "positional_bias"
    elif varid in dc_variants["low_dp"]:
        filter_col = "strand_bias"
    else:
        filter_col = "low_dp"


    # Generate info field 
    dp = cols[4]
    ad = str(int(dp) - int(cols[6])) + "," + cols[6]
    af = str(float(cols[7])/100)
    somatic_status = cols[11]

    info_col = ";".join([f'VARSTATUS={somatic_status}'])

    # Generate GT field 
    gt = "./."
    gt_col = ":".join([gt, dp, ad, af])
    out.write("\t".join([chrom, pos, ".", ref, alt, filter_col, info_col ,format_col, gt_col]) + "\n")
out.close()
