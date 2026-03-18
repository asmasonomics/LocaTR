## python 204_ltrstruc_batch_pos_extract.py [-h] ref_genome

import argparse
parser = argparse.ArgumentParser(
	description="ltrstruc_batch_pos_extract extracts positions from a reference fasta file using putative LTR retrotransposons identified by LTR_STRUC.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ref_genome", help="Fasta file, often a reference genome, searched by LTR_STRUC for LTR retrotransposons")
usr_args = parser.parse_args()

import os
import sys
import subprocess
import re
sys.path.append("full_path_to_LocaTR/00_pipeline")
import static_functions


## Identifies all files in current directory that end with ".fasta" and adds them to the seq_files list
seq_files = []
for filename in os.listdir("."):
	if filename.endswith(".fasta"):
		seq_files += [filename]

## Check for existance of BLAST database files for given reference genome. Create if they do not exist.
if os.path.isfile((usr_args.ref_genome).replace(".fa","") + ".nhr") and os.path.isfile((usr_args.ref_genome).replace(".fa","") + ".nin") and os.path.isfile((usr_args.ref_genome).replace(".fa","") + ".nsq"):
	print("\nAll BLAST database files found. Continuing to extract positions.")	
else:
	print("\nBLAST database files not present. Implementing makeblastdb protocol.")
	subprocess.call(("makeblastdb -in " + usr_args.ref_genome + " -dbtype nucl"), shell=True)
	print("BLAST database files created. Continuing to extract positions.")

final_output = open("ls_pos.txt", "w")
E = 0   #Error file counter
for filename in seq_files:
	seq_len = ((filename.split("_"))[-1]).replace(".fasta", "").replace("rc", "")   #extracts length from query file name
	seq_contig = (filename.split("_"))[0]
	if seq_contig.startswith("ch"):
			seq_contig = seq_contig.replace("ch", "")   #extracts contig/chromosome name/number and removes potential leading "ch"
	if seq_contig.startswith("0"):
			seq_contig = seq_contig.replace("0", "")   #extracts contig/chromosome name/number and removes potential leading "0"
	seq_contig = seq_contig.replace("-rc", "")

        #extracts the LTR lengths from the file headers
	header = (((open(filename).read().replace("\r", "").split("\n"))[0]).replace(">", "")).replace("-", "\t")
	
        #calls BLASTn and produces an output file "results.out" - this will be over-written on each iteration unless it is saved as an error report (below)
	subprocess.call("blastn -db " + ((usr_args.ref_genome).replace(".fa","")) + " -query " + filename + " -num_threads 8 -max_target_seqs 5 -out results.out -outfmt \"6 sacc sstart send length\"", shell=True)
	top_hit = ((open("results.out").read().split("\n"))[0]).split("\t")  #takes top hit and creates a list of the columns
        
	#if the top hit contig/chromosome is the same as the query stated AND the lengths are the same we write to file (depending on strand)
	if ((str(top_hit[0]) == str(seq_contig)) or (seq_contig == "nonchro")) and (int(top_hit[3]) == int(seq_len)):
		if (top_hit[1] < top_hit[2]):
                        #header lengths in bp are added as columns (in the else statement these are reversed like with the start/end pos so column 6 is always the most 5' LTR length - not always the 5' LTR, as -ve strand will have the 3' first)
			final_output.write(str(top_hit[0]) + "\t" + str(top_hit[1]) + "\t" + str(top_hit[2]) + "\t+\tLS\t" + header + "\n")
		else:
			final_output.write(str(top_hit[0]) + "\t" + str(top_hit[2]) + "\t" + str(top_hit[1]) + "\t-\tLS\t" + header + "\n")
        #else an error file is created suggesting to the user they should check the results file here - see README
	else:
		E += 1
		error = open("wrong_match.txt", "a")
		error.write("Error " + str(E).zfill(3) + ": " + filename + "  sequence does not return itself as the top blastn hit. Check error file.\n")
		error.close()
	
final_output.close()

pos_list = static_functions.list_initial_formatter("ls_pos.txt")
x = pos_list
x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))
i = 0
sort = None

#takes the list sorted with the lambda construct and merges positions
while i < len(x) - 1:
	# performs a fresh sort at the start of iteration if the list order has changed
	if sort == True:
		x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))
	
	if (x[i][0] == x[i + 1][0]) and (x[i][1] == x[i + 1][1]) and (x[i][2] == x[i + 1][2]):
		new = ((str(x[i])).replace("\'", "").replace("[", "").replace("]", "").replace("+", "m").replace("-", "m")).split(", ")
		new_pos = []
		count = 0
		for value in new:
			count += 1
			if count in [2, 3, 6, 7]:
				y = int(value)
				new_pos += [y]
			#not changing to int for columns which may use non numeric characters
			else:
				new_pos += [value]
                x += [new_pos]
                x.pop(i+1)
                x.pop(i) 
		sort=True 
        else:
		i += 1
		sort = False 

outfile = open("ltrstruc_positions.txt", "w")
for pos in x:
	outfile.write(str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + str(pos[3]) + "\t" + str(pos[4]) + "\t" + str(pos[5]) + "\t" + str(pos[6]) + "\n")
outfile.close()

os.remove("results.out")
os.remove("ls_pos.txt")

subprocess.call("cut -f1-5 ltrstruc_positions.txt > ls.tmp", shell=True)
subprocess.call("python full_path_to_LocaTR/001_seq_extract.py -u --prefix ltrstruc ls.tmp " + usr_args.ref_genome, shell=True)
os.remove("ls.tmp")