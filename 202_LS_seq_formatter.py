### python 202_LS_seq_formatter.py [-h] input_seq.fasta

import argparse
parser = argparse.ArgumentParser(
	description="202_LS_seq_formatter splits a given multi sequence fastA file and also reverse complements each sequence.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("input_seq", help="Fasta file to be analysed")
usr_args = parser.parse_args()

import subprocess
import sys
sys.path.append("/groups2/avian_genomes/chicken_erv/04_avian_lineage_LTR_retrotransposon_full_analysis/00_pipeline")
import static_functions

fasta_format_checker = int(subprocess.check_output("grep \"^>\" " + usr_args.input_seq + " | wc -l | awk \'{print $1}\'", shell=True))
if (fasta_format_checker == 0):
        sys.exit("Input file not in fasta format. Ensure sequences have headers starting with \">\".")
        
# Open user input file, 
seq_file = open(usr_args.input_seq).read().rstrip("\n").split(">")
seq_file.pop(0)

# extract headers and sequences, filtering out those contigs less than 2kb in length
headers = []
seq_list = []
for seq in seq_file:
	if (len(seq.split("\n")[1])>1999):
		headers += [seq.split("\n")[0]]
		seq_list += [seq.split("\n")[1]]

# split fasta into multiple files
static_functions.fast_fasta_splitter(headers, seq_list, "Y")
subprocess.call("for file in *.fas; do mv \"$file\" \"${file%.fas}.txt\"; done", shell=True)

# Run reverse complement 
i=0
for element in seq_list:
	rc = ((element.replace("A" , "t").replace("T" , "a").replace("G" , "c").replace("C" , "g")).upper())[::-1]
	out_file = open((headers[i] + "_rc.txt"), "w")
	out_file.write(">" + (headers[i]) + "_rc\n" + rc + "\n")
	out_file.close
	i+=1