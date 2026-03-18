## python 406_convert_back_to_original_contig_names.py [-h] [-p PREFIX] file_to_convert species_name_information.txt

import argparse

parser = argparse.ArgumentParser(
	description="Python script to revert position contig names back to the original in any file type (positions or sequence)",
	epilog="Author: Andrew Mason; Release: 15/06/16; Contact: andrew.mason@roslin.ed.ac.uk"
	)	
parser.add_argument("to_convert", help="positions/seq file to convert")
parser.add_argument("ref_info", help="reference genome information file from initial processing")
parser.add_argument("-p", "--prefix", default="converted", help="specify a string prefix (default = converted)")

usr_args = parser.parse_args()

# take information file and create a list of original and processed contig names
info_details = []
for line in open(usr_args.ref_info).read().rstrip("\n").split("\n"):
	info_details += [[(line.split("\t")[0]),((line.split("\t")[1]).split("_")[0]).replace("-", "_")]]

# replace contig names in file
data = open(usr_args.to_convert).read().rstrip("\n")
for pair in info_details:
	p = data.replace(pair[0], pair[1])
	data = p

ext = (usr_args.to_convert).split(".")[-1]
outfile = open(usr_args.prefix + "_original_contig_names." + ext, "w")
outfile.write(data)
outfile.close()