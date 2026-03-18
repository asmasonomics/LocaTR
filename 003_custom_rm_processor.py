## python 003_custom_rm_processor.py [-h] rm_cust_pos.txt rm_cust_seq.fasta

import argparse
import os
import sys
import subprocess

parser = argparse.ArgumentParser(
	description="custom_rm_processor takes the pos file and extracted seq file from a custom RM run and edits seq names to work with RM, performs fragment remover protocol and then sorts the new positions, and then restores the original positions.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("pos_file", help="Positions file generated from the extract_RM_positions.sh script")
parser.add_argument("seq_file", help="Sequences file generated from the above positions file")
usr_args = parser.parse_args()

# Identify the names and order of contig/chromo names and assign new, short name
switch_list = []
for line in ((subprocess.check_output(("cut -f1 " + usr_args.pos_file + " | uniq | awk \'BEGIN {i=1};{print $1\"ChR,seq\"i\"-\"; i++}\'"), shell=True)).rstrip("\n").split("\n")):
	switch_list += [[line.split(",")[0], line.split(",")[1]]]

# Read in the pos and seq files
pos_data = (subprocess.check_output(("awk \'{print $1\"ChR\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5}\' " + usr_args.pos_file), shell=True)).rstrip("\n")
seq_data = (subprocess.check_output(("cat " + usr_args.seq_file + " | rev | sed \'s/_/_RhC/3\' | rev"), shell=True)).rstrip("\n")

# Switch all the contig/chromo names for the short names assigned above
for i, j in switch_list:
	pos_data = pos_data.replace(i, j)
	seq_data = seq_data.replace(i, j)

# Write this edited data to file
pos_edit_out = ((usr_args.pos_file).split("."))[0] + "_edit.txt"
pos_edit = open(pos_edit_out, "w")
pos_edit.write(pos_data)
pos_edit.close()
seq_edit_out = ((usr_args.seq_file).split("."))[0] + "_edit.fasta"
seq_edit = open(seq_edit_out, "w")
seq_edit.write(seq_data)
seq_edit.close()
#print((sys.argv[1]).split("_")[0] + " - replacements complete")

# Perform secondary RepeatMasker with -species vertebrates -nolow flags
subprocess.call(("repeatmasker -species vertebrates -nolow " + seq_edit_out + " &> /dev/null"), shell=True)
for outfile in os.listdir("."):
	for ext in (".cat", ".cat.gz", ".masked", ".tbl"):
		try:
			os.remove(seq_edit_out + ext)
		except OSError:
			pass

# Call the fragment_remover.py protocol
subprocess.call(("python full_path_to_LocaTR/004_rm_fragment_remover.py -p " + seq_edit_out + ".out " + pos_edit_out + " &> /dev/null"), shell=True)

# Convert contig/chromo back to original
frag_data = open(((usr_args.seq_file).split("_")[0]) + "_original_pos_no_non-ltr_annotations.txt").read().rstrip("\n")
rm_data = open(seq_edit_out + ".out").read().rstrip("\n")
for i, j in switch_list:
	frag_data = frag_data.replace(j, i)
	rm_data = rm_data.replace(j, i)

# Write reverted data back to file
val_out = (usr_args.seq_file).replace("_seq.fasta", "_validated_positions.txt")
val_pos = open(val_out, "w")
val_pos.write(frag_data.replace("ChR", ""))
val_pos.close()
subprocess.call(("awk \'{if ($4!~/+|-/) {print $1\"_\"$2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6} else {print $0}}\' " + val_out + " > " + val_out.replace(".txt", ".tmp")), shell=True)
subprocess.call("sort -k1,1 -k2,2n " + val_out.replace(".txt", ".tmp") + " > " + val_out, shell=True)
subprocess.call(("python full_path_to_LocaTR/003_pos_merger.py --prefix " + ((usr_args.seq_file).split("_")[0]) + " " + val_out) + " &> /dev/null", shell=True)
rm_out = open(usr_args.seq_file + ".out", "w")
rm_out.write(rm_data.replace("ChR", ""))
rm_out.close()

# Tidy up files!
for outfile in os.listdir("."):
	if outfile.startswith((sys.argv[1]).split("_")[0]):
		for ext in ("edit.txt", "edit.fasta", "edit.fasta.out", "annotations.txt", ".tmp"):
			if outfile.endswith(ext):
				os.remove(outfile)