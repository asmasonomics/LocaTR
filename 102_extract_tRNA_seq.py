## python extract_tRNA_seq.py [-h] [--threads THREADS] <ref_genome> <tRNA_pos_file>

import argparse
parser = argparse.ArgumentParser(
	description="extract_tRNA_seq takes the output of tRNAscan-SE, groups by tRNA type (excluding \"pseudo\" annotations as these could be any tRNA), extracts sequences, performs MSA and creates profile HMMs ready for use in HMMER analysis.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ref_genome", help="Specify reference genome file searched in tRNAscan-SE")
parser.add_argument("tRNA_pos_file", help="Specify the tRNAscan-SE output file (produced by \"-o\" option)")
parser.add_argument("--threads", type=int, default=8, help="Specify number of threads (default = 8)")
usr_args = parser.parse_args()

import os
import subprocess
import sys
sys.path.append("full_path_to_LocaTR/00_pipeline")
import static_functions

startdate = subprocess.check_output("date | awk \'{print $3 \"th \" $2 \" \" $6 \"  \"  $4}\'", shell=True)

# create reference genome dict for seq extract
print("\nCreating reference genome dict for sequence extraction.")
ref_dict = static_functions.seq_dict_creator(usr_args.ref_genome)

# open tRNA output pos file and format into a list of positions
print("Formatting tRNAscan-SE output to usable positions format.")
sorted_tRNA_pos_list = subprocess.check_output("tail -n +4 " + usr_args.tRNA_pos_file + " | awk \'{print $1 \"\t\" $3 \"\t\" $4 \"\t\" $5}\' | awk \'{if($2 > $3){print $1 \"\t\" $3 \"\t\" $2 \"\t-\t\" $4} else {print $1 \"\t\" $2 \"\t\" $3 \"\t+\t\" $4}}\' | egrep -v \"Pseudo|Undet\" | sort -k5", shell=True)
print("Raw formatting complete. Processing positions into iterable lists based on coded Amino Acid Residue.")
tRNA_pos_list = static_functions.list_formatter(sorted_tRNA_pos_list, 4)

tRNA_seq_files = []
for tRNA_type in tRNA_pos_list:
        filename = (tRNA_type[0][4])[:3] + "_tRNA_seq.fasta"
        tRNA_seq_files += [filename]
        print("Extracting " + filename + ".")
        static_functions.seq_extract(filename, tRNA_type, ref_dict, "Y")

for seq_file in tRNA_seq_files:
        root_file_name = seq_file.replace("_tRNA_seq.fasta", "") 
        print("Creating MSA for " + root_file_name + " using Muscle aligner.")
        subprocess.call("muscle -in " + seq_file + " -out " + root_file_name + "_muscle_aligned.fasta > /dev/null", shell=True)
        print("Creating pHMM for aligned " + root_file_name + " sequences using HMMER hmmbuild.")
        subprocess.call("hmmbuild -o /dev/null --dna --cpu " + str(usr_args.threads) + " " + root_file_name + ".hmm " + root_file_name + "_muscle_aligned.fasta", shell=True)

os.mkdir("./02_SIE/tRNA_pHMMs")
subprocess.call("rm *_tRNA_seq.fasta *_muscle_aligned.fasta", shell=True)
subprocess.call("mv *.hmm ./02_SIE/tRNA_pHMMs/", shell=True)

end = subprocess.check_output("date | awk \'{print $3 \"th \" $2 \" \" $6 \"  \"  $4}\'", shell=True)
print("\n\nProcess complete.\nProcess began: " + str(startdate) + "Process finished: " + str(end) + "\n")