import pysam 
import argparse

parser = argparse.ArgumentParser(description='Script for testing command line arguments')
parser.add_argument('-i','--input', dest='input', help='Input delly vcf')
parser.add_argument('-f','--filter', action="store_true", dest='filter', default=False, help='Keep only variants that passed the filter')
parser.add_argument('-o','--output', dest='output', help='Output')
args = parser.parse_args()

out = open(args.output, "w")

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
            end1 = cols[pos_idx]
            if svtype == "BND":
                chrom2 = info_dc["CHR2"]
            else:
                chrom2 = chrom1 
            if svtype == "BND":
                start2 = str(int(info_dc["POS2"]) - 1)
                end2 = info_dc["POS2"]
            else:
                start2 = str(int(info_dc["END"]) - 1)
                end2 = info_dc["END"]
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
                 
            if args.filter == True:
                if filter_status == "PASS":
                    out.write("\t".join([chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, svtype, info]) + "\n")
                else:
                    pass
            else:
                out.write("\t".join([chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, svtype, info]) + "\n")

out.close()

