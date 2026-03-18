#!/bin/sh

## 305_merge_homology_positions.sh diretory/of/homology/positions ref_genome.fa
## process the output of the homology identification positions files

# help message
if [ $# -ne 2 ] || [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./305_merge_homology_positions.sh path/to/dir/with/homology_positions genome.fa"; echo ""
   echo "Please rerun script with correct usage. The script will now exit."; echo ""
   exit
fi

for i in `ls $1 | grep "positions.txt"`; do echo -e "\n$i"; tot=`awk '{sum+=(($3+1)-$2)}END{print sum}' $i`; printf "Total: %'d bp\n" $tot; done; echo ""

echo "Merging positions files"
cat *positions.txt | sort -k1,1 -k2,2n | uniq > homology_ID_concat.tmp
python full_path_to_LocaTR/002_pos_merger.py --prefix homology_ID homology_ID_concat.tmp &> /dev/null
rm homology_ID_concat.tmp

echo -e "\nExtracting sequences\n"
python full_path_to_LocaTR/001_seq_extract.py -u --prefix homology_ID homology_ID_merged_positions.txt $2

tot=`awk '{sum+=(($3+1)-$2)}END{print sum}' homology_ID_merged_positions.txt`
printf "Merged total: %'d bp\n" $tot; echo ""