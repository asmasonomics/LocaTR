## python 201_extract_LH_positions.py [-h] species_name_LH_out.sorted.gff3 species_name_process.fa

import argparse
import subprocess

parser = argparse.ArgumentParser(
	description="extract_LH_positions.py takes the output from LTR Harvest (sorted gff3 file) and extracts the positions and sequences",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("LH_out", help="sorted gff3 file")
parser.add_argument("ref", help="reference genome for sequence extraction")
usr_args = parser.parse_args()

# get length of zero indexed string name for altering naming below
zero_length = len(((subprocess.check_output("grep -m 1 \">\" " + usr_args.ref, shell=True)).rstrip("\n")).replace(">seq",""))

# remove all lines starting with # and extract only those detailing the LTR retrotransposon position. Then return the sequence name to that before the suffix array
subprocess.call("grep -v \"^#\" " + usr_args.LH_out + " | awk '$3==\"LTR_retrotransposon\"{print $1 \"\t\" $4 \"\t\" $5}' | awk '{var=\"%0" + str(zero_length) + "d\"}{print \"seq\"(sprintf(var,(substr($1,4,length($1))+1))) \"\t\" $2 \"\t\" $3 \"\t+\tLH\"}' | sort -k1,1 -k2,2n > " + (usr_args.LH_out).split(".")[0] + "_positions.txt", shell=True)

# extract sequences
subprocess.call("python full_path_to_LocaTR/001_seq_extract.py -u --prefix " + (usr_args.LH_out).split(".")[0] + " " + (usr_args.LH_out).split(".")[0] + "_positions.txt " + usr_args.ref, shell=True)