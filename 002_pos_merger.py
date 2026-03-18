## python 002_pos_merger.py [-h] [--prefix PREFIX] pos_file 

import argparse
parser = argparse.ArgumentParser(
	description="pos_merger takes a file of positions, sorts them and then merges.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("pos_file", help="Positions file to be analysed - column order (tab spaced): chr st end str source_program")
parser.add_argument("--prefix", help="Specify your own output file prefix (default = None)")
usr_args = parser.parse_args()

import sys
import subprocess
sys.path.append("full_path_to_LocaTR/00_pipeline")
import static_functions

x = static_functions.list_initial_formatter(usr_args.pos_file)
tmp_log = open("log.tmp", "w")
merged_list = static_functions.positions_merger(x, tmp_log)
try:
        subprocess.call("rm log.tmp", shell=True)
except OSError:
        pass

outname = "merged_positions.txt"
if usr_args.prefix:
        outname = usr_args.prefix + "_merged_positions.txt"
outfile = open(outname, "w")
for pos in merged_list:
	outfile.write(str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + str(pos[3]) + "\t" + str(pos[4]) +"\n")
outfile.close()