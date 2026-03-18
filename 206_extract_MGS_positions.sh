# Usage: ./206_extract_MGS_positions.sh <MGS_ltr.out> <ref_genome.fa>

## process the MGEScan_LTR output file to a more usable format
# Script reads file and excludes blank lines and 'clade lines'
# If 'start pos' values are less than 1 these are replaced with '1'
# Chromosome names are formatted to not have their file extension (will deal with .fa*)
# File output in standard positions ordering

#!/bin/sh

# help message
if [ $# -ne 2 ] || [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./206_extract_MGS_positions.sh ltr.out genome.fa"; echo ""
   echo "Please rerun script with correct usage. The script will now exit."; echo ""
   exit
fi

# set output prefix
outname=`echo $2 | rev | cut -d/ -f1 | rev | cut -d. -f1`

# from ltr.out file, exclude blank lines and "clade" lines, sort and process lines where 
# the start value is less than 1. Sorts contig names, process output file and sorts.
grep -v "^$" $1 \
| grep -v '^[0-9][0-9]*--' \
| sort -n \
| awk '{OFS="\t"}$2<1{$2="1"}1' \
| awk '{OFS="\t"}{end=index($1,".fa")}{$1=substr($1,0,end-1)}1' \
| awk '{print $1 "\t" $2 "\t" $5 "\t" $6 "\tMGS"}' \
| sort -k1,1 -k2,2n > "${outname}_MGS_positions.txt"

# extract sequences
python full_path_to_LocaTR/001_seq_extract.py -u --prefix "${outname}_MGS" "${outname}_MGS_positions.txt" $2