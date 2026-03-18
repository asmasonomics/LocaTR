#!/bin/sh

## ./404_merge_positions.sh all_files_to_merge
## takes a list of all positions files to merge, concatenates and sorts them, and then runs the pos_merger protocol.

if [ $1 == "-h" ] || [ $1 == "-help" ]; then
   echo ""; echo "Usage: ./404_merge_positions.sh all_files_to_merge"; echo ""
   echo "Please rerun script with correct usage. The script will now exit."; echo ""
   exit
fi

cat $* | sort -k1,1n -k2,2n | uniq > concat_pos.tmp

python full_path_to_LocaTR/002_pos_merger.py --prefix all concat_pos.tmp

rm concat_pos.tmp