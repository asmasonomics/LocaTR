## python 405_secondary_BLAST_analysis.py [-h] [-p PREFIX] SIE_positions.txt Homology_positions.txt /full/path/to/ref_genome.fa

import argparse
import os
import subprocess


parser = argparse.ArgumentParser(
	description="For completing the secondary BLAST analysis of the reference genome following SIE and homology identification, validation  and position list merging",
	epilog="Author: Andrew Mason; Release: 15/06/16; Contact: andrew.mason@roslin.ed.ac.uk"
	)	
parser.add_argument("sie_pos", help="positions file from SIE merge")
parser.add_argument("hom_pos", help="positions file from Homology merge")
parser.add_argument("ref_gen", help="reference genome with full path, no relative")
parser.add_argument("-p", "--prefix", default="specify a string prefix", help=" (default = USER)")

usr_args = parser.parse_args()


# Part 1 - Overlap SIE positions list with Homology list and keep only those SIE positions not found within the Homology list

print("\nReformatting positions files as BED files")
subprocess.call("awk \'{print $1 \"\t\" ($2-1) \"\t\" $3 \"\t\" $5 \"\t0\t\" $4}\' " + usr_args.sie_pos + " > " + (usr_args.sie_pos).replace(".txt", ".bed"), shell=True)
subprocess.call("awk \'{print $1 \"\t\" ($2-1) \"\t\" $3 \"\t\" $5 \"\t0\t\" $4}\' " + usr_args.hom_pos + " > " + (usr_args.hom_pos).replace(".txt", ".bed"), shell=True)
print("\nOverlapping SIE and Homology positions list, keeping only those unique to SIE list.")
subprocess.call("intersectBed -v -a " + (usr_args.sie_pos).replace(".txt", ".bed") + " -b " + (usr_args.hom_pos).replace(".txt", ".bed") + " | awk \'{print $1 \"\t\" ($2+1) \"\t\" $3 \"\t\" $6 \"\tuSIE\"}\' > " + usr_args.prefix + "_unique_SIE_positions.txt", shell=True)
os.remove((usr_args.sie_pos).replace(".txt", ".bed"))
os.remove((usr_args.hom_pos).replace(".txt", ".bed"))
print("Extracting unique SIE sequences.")
subprocess.call("python full_path_to_LocaTR/001_seq_extract.py --prefix " + usr_args.prefix + " -u " + usr_args.prefix + "_unique_SIE_positions.txt " + usr_args.ref_gen, shell=True)

# Part 2 - use the extracted sequences as the new reference list for a refBLAST rerun

print("\nPerforming secondary refBLASTsearch style analysis using unique SIEs as reference sequences")
cwd = (os.getcwd()) + "/"
subprocess.call("python full_path_to_LocaTR/302_refBLASTsearch.py " + cwd + usr_args.prefix + "_seq.fasta " + usr_args.ref_gen + " &> sec_refBLAST.log", shell=True)
os.remove(usr_args.prefix + "_seq.fasta")
os.remove(usr_args.prefix + "_unique_SIE_positions.txt")
subprocess.call("mv refBLAST_final_identified_sequences.fasta " + usr_args.prefix + "_sec_refBLAST_sequences.fasta", shell=True)
subprocess.call("mv refBLAST_final_identified_positions.txt " + usr_args.prefix + "_sec_refBLAST_positions.txt", shell=True)
subprocess.call("mv refBLAST_files sec_refBLAST_files", shell=True)
print("Done.\n")