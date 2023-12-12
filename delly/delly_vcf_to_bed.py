import argparse
import os

parser = argparse.ArgumentParser(description='Convert VCF from delly to browser extendable format')
parser.add_argument('-i','--input', dest='input', help='Input delly vcf')
parser.add_argument('-s','--samplename', dest='samplename', help='Samplename')
parser.add_argument('-o','--outdir', dest='outdir', help='Outdir')
args = parser.parse_args()

'''
Script converts the vcf-file into bed and bedpe formats : 
1. Deletions => bed 
2. Inversion => bed  
3. Tandem duplication => bed
4. Translocation => bedpe
'''

# Output 
out_del = open(os.path.join(args.outdir, args.samplename + "_deletions.bed"), "w")
out_inv = open(os.path.join(args.outdir, args.samplename + "_inversions.bed"), "w")
out_dup = open(os.path.join(args.outdir, args.samplename + "_dup.bed"), "w")
out_bnd = open(os.path.join(args.outdir, args.samplename + "_translocations.bedpe"), "w")

# Iterate over input lines 
with open(args.input) as infile:
    for ln in infile:
        if ln.startswith("#"):
            if ln.startswith("#CHROM"):
                cols = ln.split("\t")
                chrom_idx = cols.index("#CHROM")
                pos_idx = cols.index("POS")
                id_idx = cols.index("ID")
                ref_idx = cols.index("REF")
                alt_idx = cols.index("ALT")
                qual_idx = cols.index("QUAL")
                filter_idx = cols.index("FILTER")
                info_idx = cols.index("INFO")
        else:
            # Split col 
            cols = ln.split("\t")
                
            # Find if sv site location is PRECISE or not 
            sv_site_location = cols[info_idx].split(";")[0]
            precise = "Yes" if sv_site_location == "PRECISE" else "No"

            # Parse info field as dictionary
            info_dc = dict([tuple(item.split("=")) for item in cols[info_idx].split(";")[1:]])
            
            # Get SV type 
            svtype = info_dc["SVTYPE"]
            
            # Construct the BEDPE fields
            chrom1 = cols[chrom_idx]
            start1 = str(int(cols[pos_idx]) - 1)

            if svtype == "BND":

                # End 
                end1 = cols[pos_idx]
                
                # Chromosome 2 
                chrom2 = info_dc["CHR2"]

                # Start 2
                start2 = str(int(info_dc["POS2"]) - 1)

                # End 2
                end2 = info_dc["POS2"]

            else:
                
                end1 = info_dc["END"]

            name = cols[id_idx]
            score = "."
            strand1 = "."
            strand2 = "."
                
            # Custom info column
            filter_status = cols[filter_idx]

            # Extract mapping qualities
            mapq = info_dc["MAPQ"]

            # Paired-end signature induced connection type
            CT = info_dc["CT"]

            if "PE" in info_dc.keys():
                pe = info_dc["PE"]
            else:
                pe = "0"
                    
            if "SR" in info_dc.keys():
                sr = info_dc["SR"]
            else:
                sr = "0"
                
            if "SRMAPQ" in info_dc.keys():
                srmapq = info_dc["SRMAPQ"]
            else:
                srmapq = "NA"

            # Info string 
            info = f'FILTER={filter_status},PRECISE={precise},MAPQ={mapq},SRMAPQ={srmapq},CT={CT},PE={pe},SR={sr}'

            if svtype == "BND":
                # IF SV is a translocation produce bedpe
                out_bnd.write("\t".join([chrom1, start1, end1, 
                                                 chrom2, start2, end2, 
                                                 name, score, strand1, 
                                                 strand2, svtype, info]) + "\n")
            elif svtype == "DEL":
                # If SV is a deletion produce bed 
                out_del.write("\t".join([chrom1, start1, end1, 
                                                 name, score, strand1, 
                                                 svtype, info]) + "\n")
                        
            elif svtype == "DUP":
                # If SV is duplication produce bed 
                out_dup.write("\t".join([chrom1, start1, end1, 
                                                 name, score, strand1, 
                                                 svtype, info]) + "\n")
            else:
                # If SV is an inversion produce bed
                out_inv.write("\t".join([chrom1, start1, end1, 
                                                 name, score, strand1, 
                                                 svtype, info]) + "\n")
out_del.close()
out_inv.close()
out_dup.close()
out_bnd.close()

