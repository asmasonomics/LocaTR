## python 304_extract_dirs_positions.py [-h] [-p PREFIX] annotCoGenom.tmp

import argparse
import subprocess

parser = argparse.ArgumentParser(
        description="extract_dirs_positions takes the annotCoGenom file from ReDoSt and extracts positions.",
        epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("redost", help="annotCoGenom.tmp file")
parser.add_argument("-p", "--prefix", default="species", help="prefix for outfile, usually species name (default = species)")
usr_args = parser.parse_args()

# create output file
dirs_out = open(usr_args.prefix + "_dirs_positions.txt","w")

# for each line in the DIRS output, extract the domain e values as floats
for line in open(usr_args.redost).read().rstrip("\n").split("\n"):
	
	rt_eval = float(line.split(" ")[3])
	mt_eval = float(line.split(" ")[6])
	yr_eval = float(line.split(" ")[-1])
	
	# if the RT e value is less than 1e-15 and either MT or YR e-values are less than 1e-12, extract positions
	if (rt_eval < float("1e-15")) and (mt_eval < float("1e-12") or yr_eval < float("1e-12")):
		# define contig, first and last rt pos, and last yr pos
		contig = (line.split(" ")[0]).split("@")[0]
		rt_fp = int(line.split(" ")[1])
		rt_lp = int(line.split(" ")[2])
		yr_lp = int(line.split(" ")[-2])
		
		# check orientation of elements and write accordingly
		# add 350bp of flanking DNA on either side for putative inverted repeats
		if rt_fp < rt_lp:
			dirs_out.write(contig + "\t" + str(rt_fp - 350) + "\t" + str(yr_lp + 350) + "\t+\tDIRS\n")
		else:
			dirs_out.write(contig + "\t" + str(yr_lp - 350) + "\t" + str(rt_fp + 350) + "\t-\tDIRS\n")
dirs_out.close()

if int(subprocess.check_output("wc -l " + usr_args.prefix + "_dirs_positions.txt | awk \'{print $1}\'", shell=True).rstrip("\n")) == 0:
	print("No validated DIRS positions in " + usr_args.prefix + "_dirs_positions.txt")
	subprocess.call("rm " + usr_args.prefix + "_dirs_positions.txt", shell=True)