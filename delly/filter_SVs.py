import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Filter SVs')
parser.add_argument('-i','--input', dest='input', help='Input delly vcf')
parser.add_argument('-f','--filter_pass', action="store_true", dest='filter_pass', default=False, help='Keep only variants that passed the Delly quality filtering')
parser.add_argument('-p','--filter_imprecise', action="store_true", dest='filter_imprecise', default=False, help='Keep only variants that are considered precise')
parser.add_argument('-e','--pe_support', dest='pe_support', default = 0, help='Keep variants with supporting PE reads equal or higher given threshold')
parser.add_argument('-r','--sr_support', dest='sr_support', default = 0, help='Keep variants with supporting SR reads equal or higher given threshold')
parser.add_argument('-s','--total_support', dest='total_support', default = 0, help='Keep variants with total supporting reads (SR + PE) more than given threshold')
parser.add_argument('-o','--output', dest='output', help='output')
args = parser.parse_args()

# Output 
out = open(args.output, "w")

# Iterate over lines 
with open(args.input) as infile:
    for i,ln in enumerate(infile):
        if i == 0:
            cols = ln[:-1].split("\t")
            info_col_idx = cols.index("INFO")
            total_support_idx = cols.index("TOTAL_SUPPORTING_READS")
            out.write(ln)
        else:
            cols = ln[:-1].split("\t")
            info_col = cols[info_col_idx]
            info_dc = dict([tuple(item.split("=")) for item in info_col.split(",")])

            # Get info fields to consider in filtering 
            pass_filter = info_dc["FILTER"]
            imprecise_filter = info_dc["PRECISE"]
            pe_support = int(info_dc["PE"])
            sr_support = int(info_dc["SR"])
            total_support = float(cols[total_support_idx])

            filters = []

            # Consider Delly quality filter
            if args.filter_pass == True:
                if pass_filter == "PASS":
                    filters.append("PASS")
                else:
                    filters.append("FAIL")
            else:
                filters.append("PASS")

            # Consider whether localisation of SV is 
            # precise
            if args.filter_imprecise == True:
                if imprecise_filter == "Yes":
                    filters.append("PASS")
                else:
                    filters.append("FAIL")
            else:
                filters.append("PASS")
            
            # Consider the number of PE reads supporting 
            # SV
            if pe_support >= int(args.pe_support):
                filters.append("PASS")
            else:
                filters.append("FAIL")

            # Consider the total number of reads supporting 
            # SV
            if sr_support >= int(args.sr_support):
                filters.append("PASS")
            else:
                filters.append("FAIL")
            
            # Consider the total number reads supporting 
            # SV (PE + SR)
            if total_support >= int(args.total_support):
                filters.append("PASS")
            else:
                filters.append("FAIL")

            # If any of the filters fail do not report SV 
            if any(np.array(filters) == "FAIL"):
                pass 
            else:
                out.write(ln)
out.close()
