library(argparse)
library(StructuralVariantAnnotation)

# Command line params 
parser <- ArgumentParser()

# by default ArgumentParser will add an help option
parser$add_argument("-v", "--vcf", help="Input VCF")
parser$add_argument("-f", "--filter", action="store_true", default = F, help="Keep only passed variants")
parser$add_argument("-o", "--output", help="Output bedpe file")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# Read in VCF
vcf <- VariantAnnotation::readVcf(args$vcf)

# Convert to breakpoint ranges format
break.point.ranges <- breakpointRanges(vcf)

if (args$filter == T){
   # Only keep SVs with FILTER == PASS 
   break.point.ranges = break.point.ranges[break.point.ranges$FILTER == "PASS",]
}

# Conversion to bedpe-format 
rtracklayer::export(breakpointgr2pairs(break.point.ranges), con=args$output)
