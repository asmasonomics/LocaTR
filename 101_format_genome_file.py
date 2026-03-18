## python 101_format_genome_file.py [-h] genome.fa

from __future__ import division
import argparse
import re
import subprocess

parser = argparse.ArgumentParser(
	description="format_genome_file.py takes the specified reference genome file and creates a formatted version ready for analysis",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ref", help="genome file in fasta format")
usr_args = parser.parse_args()

###################
#
# Performance:
# Takes less than 1 minute to process galGal4 assembly. File is 1017Mb, containing 1.073 Gbp of sequence over 15932 contigs.
#
###################

# unzip the genome if it is a zipped file. NOTE - this does remove the zipped file
if (usr_args.ref).endswith(".gz"):
	subprocess.call("gunzip -c " + usr_args.ref + " > " + (usr_args.ref).replace(".gz",""), shell=True)

# read in the genome, separate headers and sequence, get sequence on single lines
genome = open((usr_args.ref).replace(".gz","")).read().rstrip("\n")
# if the genome file was zipped, remove the unzipped version
if (usr_args.ref).endswith(".gz"):
	subprocess.call("rm " + (usr_args.ref).replace(".gz",""), shell=True)

headers = []
for line in genome.split("\n"):
    if line.startswith(">"):
        headers += [(re.sub(r"[ ]+", "_", (line.replace(">", "").replace("_","-"))))]
header_count = len(headers)

sequences = re.split(r"_", (re.sub(r"[\r|\n]", "", (re.sub(r"(>(.+)\n)", "_", genome)))).upper())    
sequences.pop(0)

# create new genome file and reference text document
i=1
prefix = ((usr_args.ref).split("/")[-1]).split(".")[0]
genome_out = open(prefix + "_processed.fa", "w")
genome_ref = open(prefix + "_information.txt", "w")

for header in headers:
    new_header = ">seq" + str((str(i).zfill(len(str(header_count)))))
    genome_out.write(new_header + "\n" + (sequences[i-1]) + "\n")

    seq_length = len(sequences[i-1])
    N_count = (sequences[i-1]).count("N")
    seq_GC = round((float(100/(seq_length - N_count)) * ((sequences[i-1]).count("G") + (sequences[i-1]).count("C"))), 2)
    seq_Ns = round((float(100/(seq_length)) * (N_count)), 2)
    genome_ref.write(new_header.replace(">", "") + "\t" + headers[i-1] + "\t" + str(seq_length) + "\t" + str(seq_Ns) + "\t" + str(seq_GC) + "\n")

    i+=1

genome_out.close()
genome_ref.close()

# create limited directory structure
subprocess.call("mkdir 01_Homology 02_SIE 03_Validated", shell=True)
