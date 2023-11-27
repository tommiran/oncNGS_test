from intervaltree import Interval, IntervalTree
import argparse

parser = argparse.ArgumentParser(description='Annotate translocation sites')
parser.add_argument('-i','--input', dest='input', help='Input bedpe')
parser.add_argument('-g','--gtf', dest='gtf', help='Input GTF')
parser.add_argument('-o','--output', dest='output', help='Output')
args = parser.parse_args()

"""
Functions 
""" 
#
# Function prepares are dictionary of genes based 
# on the gene-features extracted from the given 
# GTF file. 
#
def prepareGeneDatabase(gtf_file, target_chromosomes):

    # Dictionary of gene structures 
    gtf_dict = {}

    # Interval dictionary 
    intervals_dc = {}

    # Iterate over GTF lines
    with open(gtf_file) as infile:
        for ln in infile:
            if ln.startswith("#!"):
                pass
            else:
                cols = ln[:-1].split("\t")
        
                # Extract information from columns
                feature = cols[2]
                chrom = "chr" + cols[0]
                start = str(int(cols[3]) - 1)
                end = cols[4]

                # Only consider the regular chromsomes 
                if chrom in target_chromosomes:
        
                    if feature == "gene":
                
                        # Parse last column by the delim (",") as a
                        # dictionary
                        ls = []
                        for item in cols[-1].split(";")[:-1]:
                            key,value = item.split()[0].lstrip(),item.split()[1].replace("\"","")
                            ls.append((key, value))
                        info_dc = dict(ls)
        
                        # Make a key to gtf dc for east query
                        gtf_dc_key = chrom + ":" + start + "-" + end
        
                        # Check if already in dc
                        if gtf_dc_key not in gtf_dict.keys():
            
                            # Initiate new list of entries
                            gtf_dict[gtf_dc_key] = []
            
                            # Prepare new entry 
                            gtf_entry =  {"Feature" : feature, 
                                          "Coordinates" : chrom + ":" + start + "-" + end,
                                          "Gene_ID" : info_dc.get("gene_id"),
                                          "Gene_name" : info_dc.get("gene_name"), 
                                          "Gene_biotype" : info_dc.get("gene_biotype")}        
            
                            gtf_dict[gtf_dc_key].append(gtf_entry)
                        else:
                            # Prepare new entry 
                            gtf_entry =  {"Feature" : feature, 
                                          "Coordinates" : chrom + ":" + start + "-" + end,
                                          "Gene_ID" : info_dc.get("gene_id"),
                                          "Gene_name" : info_dc.get("gene_name"), 
                                          "Gene_biotype" : info_dc.get("gene_biotype")}        
            
                            gtf_dict[gtf_dc_key].append(gtf_entry)           
            
                        # Append interval 
                        if chrom in intervals_dc.keys():
                            intervals_dc[chrom].append((int(start), int(end)))
                        else:
                            intervals_dc[chrom] = []
                    else:
                        pass
                else:
                    pass

    # Generate an interval object 
    # with chromosomes as the keys 
    # and interval trees as values
    intervals_obj = {}

    for key in intervals_dc.keys():
        intervals_obj[key] =  IntervalTree.from_tuples(intervals_dc[key])

    # Generate database object
    # with gene information and interval trees
    gene_db_obj = {"gene_info" : gtf_dict, 
                   "intervals" : intervals_obj}

    return(gene_db_obj)


#
# Prepare annotation entry for the output 
#
def prepAnnotEntry(x):
    annot_entry = []
    for item in x:
       for key in item.keys():
            if item[key] != None:
                annot_entry.append(key + "=" + item[key])
            else:
                annot_entry.append(key + "=" + "NA")
    return(",".join(annot_entry))


# 
# Fetch the genes associated with overlapped 
# intervals 
#
def getGeneInfo(chrom, overlapped_intervals, gene_database):

    # Collect all overlapped genes
    overlapped_gene_annotations = []

    for interval in overlapped_intervals:
        qkey = chrom + ":" + str(interval[0]) + "-" + str(interval[1])
        match_genes = gene_database["gene_info"][qkey]
        for match_gene in match_genes:
            overlapped_gene_annotations.append(match_gene)

    return(overlapped_gene_annotations)


#
# Helper function for checking whether intervals overlap
#
def intervalsOverlap(interval1, interval2):
    
    # Get chromosomes  
    chrom1 = interval1.split(":")[0]
    chrom2 = interval2.split(":")[0]
    if chrom1 == chrom2:

        # As chromosomes match check positional coordinates
        start1 = int(interval1.split(":")[1].split("-")[0])
        end1 = int(interval1.split(":")[1].split("-")[1])
        start2 = int(interval2.split(":")[1].split("-")[0])
        end2 = int(interval2.split(":")[1].split("-")[1])

        if ((end2 < start1) & (start2 > end1)):
            return(False)
        else:
            return(True)
    else:
        return(False)
    
#
# Find all possible fusion partner-pairs
# These genes must not overlap in the absense 
# of sv
#
def findFusionPartners(overlapped_gene_annotations_left, overlapped_gene_annotations_right): 

    # Fusion pairs to be reported
    fusion_pairs = []

    for ov_left in overlapped_gene_annotations_left:
        # Extract the gene name and coordinates 
        coord_left = ov_left["Coordinates"]
        gene_left = ov_left["Gene_name"]

        # In the case that there is no gene symbol we 
        # use ENSEMBL id instead
        if gene_left == None:
            gene_left = ov_left["Gene_ID"]

        for ov_right in overlapped_gene_annotations_right:
            # Extract the gene name and coordinates 
            coord_right = ov_right["Coordinates"]
            gene_right = ov_right["Gene_name"]

            if gene_right == None:
                gene_right = ov_right["Gene_ID"]

            if ((coord_left != None) & (coord_right != None)):
                # Check overlap 
                if intervalsOverlap(coord_left, coord_right) == False:
                    # Intervals do not overlap so add pair to the list of 
                    # fusion-pairs 
                    fusion_pairs.append(gene_left + ":" + gene_right)
                else:
                    # Intervals overlap => ignore
                    pass
            else:
                pass 

    return(fusion_pairs)

"""
Main 
"""

# Harcoded target chromosomes 
target_chromosomes = ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"]

# Prepare gene database 
gene_database = prepareGeneDatabase(args.gtf, target_chromosomes)

# Open file for output 
out = open(args.output, "w")

# Iterate over input BEDPE file
header = ["#CHROM","START1","END1","CHROM2","START2","END2", "SCORE", "STRAND1", "STRAND2","SVTYPE","INFO","Predicted fusions"]
out.write("\t".join(header) + "\n")

with open(args.input) as infile:
    for ln in infile:
        cols = ln[:-1].split("\t")
        chrom, start, end, chrom2, start2, end2 = cols[:6] 
    
        # Test "left" interval 
        if chrom in target_chromosomes:
            results_left = gene_database["intervals"][chrom][int(start):int(end)]
        else:
            results_left = []
        
        # Test "right" interval
        if chrom2 in target_chromosomes:
            results_right = gene_database["intervals"][chrom2][int(start2):int(end2)]
        else:
            results_right = []
        
        # Fetch all the overlapped Genes
        annotation_left = getGeneInfo(chrom, results_left, gene_database)

        annotation_right = getGeneInfo(chrom2, results_right, gene_database)

        # Find fusion partners 
        # The genes must not overlap in the absense of the sv
        fusion_pairs= findFusionPartners(annotation_left, annotation_right)
    
        if fusion_pairs != []:
            out.write("\t".join(cols) + "\t" + ",".join(fusion_pairs) + "\n") 
        else:
            pass
    
out.close()
