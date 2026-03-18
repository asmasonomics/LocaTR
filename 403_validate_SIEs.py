### python 403_validate_SIEs.py [-h] [--dna_hmms DNA_HMMS] [--protein_hmms PROTEIN_HMMS] [--trna_hmms TRNA_HMMS] [--out_prefix OUT_PREFIX] [--threads THREADS] [--species SPECIES] positions_file seq_file

## defining user argument parsing

from __future__ import division
import argparse

parser = argparse.ArgumentParser(
	description="Python script for validating putative LTR retrotransposons using RepeatMasker and profile Hidden Markov Models for both nucleotide and protein sequences.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk"
	)

required_args_group = parser.add_argument_group(title="Required argument definitions", description="Required user defined arguments. Strict order is required.")
required_args_group.add_argument("positions_file", help="Specify the positions file")
required_args_group.add_argument("seq_file", help="Specify extracted sequence file")

run_options_group = parser.add_argument_group(title="Run time options", description="User specifiable run time options. Deafult values provided.")
run_options_group.add_argument("-d", "--dna_hmms", default="./401_Dfam_nucl_pHMMs/", help="Specify directory containing DNA pHMMs (default = ./401_Dfam_nucl_pHMMs/)")
run_options_group.add_argument("-p", "--protein_hmms", default="./402_GyDB_protein_pHMMs/", help="Specify directory containing protein pHMMs (default = ./402_GyDB_protein_pHMMs/)")
run_options_group.add_argument("-t", "--trna_hmms", default="./tRNA_pHMMs/", help="Specify directory containing tRNA pHMMs (default = ./tRNA_pHMMs/)")
run_options_group.add_argument("-op", "--out_prefix", default="validation" ,help="Specify output file prefix (default = validation)")
run_options_group.add_argument("--threads", type=int, default=8, help="<INT> Specify the number of threads (default = 8)")
run_options_group.add_argument("--species", default="vertebrates", help="Specify your RepeatMasker target species name - see RepBase.org for suitable names (default is \"vertebrates\").")

usr_args = parser.parse_args()

import sys
import subprocess
import os
import csv
import re
from collections import OrderedDict

startdate = subprocess.check_output("date | awk \'{print $3 \"th \" $2 \" \" $6 \"  \"  $4}\'", shell=True)

sys.path.append("full_path_to_LocaTR/00_pipeline")
import static_functions

prefix = usr_args.out_prefix
num_threads = usr_args.threads

if os.path.isfile(prefix + "_validation_log_file.txt"):                      # log_file.txt initiation
	log_file = open(prefix + "_validation_log_file.txt", "a")
	log_file.write("\n\n\n\n\n------------------------------------------------------------------------------------------\n\n\n" + str(startdate))
else:
	log_file = open(prefix + "_validation_log_file.txt", "w")
	log_file.write(str(startdate))

## Part 1 - process positions file, extract sequences

hit_num_before_validation = subprocess.check_output("wc -l " + usr_args.positions_file +" | awk \'{print $1}\'", shell=True)
total_hit_length_before_validation = subprocess.check_output("awk \'{sum+=($3-$2)}END{print sum / 1000000}\' " + usr_args.positions_file, shell=True)

# Processing sequence file
print("Processing putative LTR retrotransposon sequences and splitting into individual query files.");log_file.write("Processing putative LTR retrotransposon sequences and splitting into individual query files.")
indiv_seq_files = static_functions.fasta_splitter((open(usr_args.seq_file).read().rstrip("\n")), "N")

## Part 2 - Prepare files and conduct HMMER analysis
print("Processing HMM profiles to create a HMM database flatfile and a binary indexed HMM database.");log_file.write("Processing HMM profiles to create a HMM database flatfile and a binary indexed HMM database.")

# DNA HMMs
print("Processing nucleotide HMMs.");log_file.write("Processing nucleotide HMMs.")
static_functions.run_hmm_press(usr_args.dna_hmms, "nucl_pHMMs", log_file)

# tRNA HMMs
print("Processing tRNA HMMs.");log_file.write("Processing tRNA HMMs.")
static_functions.run_hmm_press(usr_args.trna_hmms, "tRNA_pHMMs", log_file)

# Protein HMMs
print("Processing protein HMMs.");log_file.write("Processing protein HMMs.")
static_functions.run_hmm_press(usr_args.protein_hmms, "prot_pHMMs", log_file)


# prepare output files
seq_info_file = open("seq_info.txt", "w")
seq_info_file.write("--Putative_element_info--\t-----\t-----\t-----\t-----\t-----\nChromosome\tStart\tEnd\tStrand\tLength\tID_prog\n")
rm_file = open("repeatmasker_scores.txt", "w")
rm_file.write("--RepeatMasker--\t-----\t-----\t-----\nAll_matches\tLTR_coverage\tLTR_position\tLTR_strand\n")
nucl_file = open("validated_nucl_scores.txt", "w")
nucl_file.write("--DNA_motifs--\t-----\t-----\t-----\t-----\t-----\nTotal_annotations\tLTR_related_annotations\tBest_hit\tStrand\tE-value\tBest_hit_position\n")
tRNA_file = open("validated_tRNA_scores.txt", "w")
tRNA_file.write("--tRNA--\t-----\t-----\nResidue\tE-value\tStrand\n")
prot_file = open("validated_prot_scores.txt", "w")
prot_file.write("--GAG--\t-----\t-----\t-----\t--Protease--\t-----\t-----\t-----\t--INT--\t-----\t-----\t-----\t--RT--\t-----\t-----\t-----\t--RNaseH--\t-----\t-----\t-----\t--dUTPase--\t-----\t-----\t-----\t--ENV--\t-----\t-----\t-----\t--ATF--\t-----\t-----\t-----\t--CHR--\t-----\t-----\t-----\t--MOV--\t-----\t-----\t-----\t--TAV--\t-----\t-----\t-----\t--VAP--\t-----\t-----\t-----\t--Accessory_proteins--\t-----\t-----\t-----\nMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\tMatches\tBest_match\tE-value\tStrand\n")

def prot_domain_sorter(domain, hmm_out_file, last_pos):
        report_line = ""
        if (int(subprocess.check_output("egrep \"" + domain + "\" " + hmm_out_file + " | wc -l", shell=True)) > 0):
                total_matches = (subprocess.check_output("egrep \"" + domain + "\" " + hmm_out_file + " | wc -l", shell=True)).rstrip("\n")
                if (last_pos == True):
                        best_match = (subprocess.check_output("egrep \"" + domain + "\" " + hmm_out_file + " | sort -g -k7 | head -n 1 | awk \'{print $1}\'", shell=True)).rstrip("\n")
                else:
                        best_match = ((subprocess.check_output("egrep \"" + domain + "\" " + hmm_out_file + " | sort -g -k7 | head -n 1 | awk \'{print $1}\'", shell=True)).rstrip("\n").split("_"))[1]
                e_val = (subprocess.check_output("egrep \"" + domain + "\" " + hmm_out_file + " | sort -g -k7 | head -n 1 | awk \'{print $7}\'", shell=True)).rstrip("\n")
                strand = (hmm_out_file.split("_"))[2]
                if (int(((subprocess.check_output("egrep \"" + domain + "\" " + hmm_out_file + " | sort -g -k7 | head -n 1 | awk \'{print $4}\'", shell=True)).rstrip("\n").split("_"))[-1]) > 3):
                        if (strand == "-"):
                                strand = "+"
                        else:
                                strand = "-"
                report_line = (str(total_matches) + "\t" + best_match + "\t" + str(e_val) + "\t" + strand + "\t")    
        else:
                report_line = "0\tNULL\tNULL\tNULL\t"
                
        return report_line
        

for seq_file in indiv_seq_files:
        ## Produce putative sequence info file
        details = seq_file
        underscore_count = ((seq_file.split("-"))[0]).count("_")
        wanted_underscores = int(underscore_count) - 1
        if (wanted_underscores > 0):
                details = re.sub(r"_", r"-", seq_file, wanted_underscores)
        chromo = (str((details.split("_"))[0])).replace("-", "_")
        start_pos = str((((details.split("_"))[1]).split("-"))[0])
        end_pos = str((((details.split("_"))[1]).split("-"))[1])
        strand = str((details.split("_"))[2])
        ID_prog = ((str((details.split("_"))[3])).split("."))[0]
        seq_length = (int(end_pos) + 1) - int(start_pos)
        seq_info_file.write(chromo + "\t" + start_pos + "\t" + end_pos + "\t" + strand + "\t" + str(seq_length) + "\t" + ID_prog + "\n")
        


        ## Nucleotide HMM analysis
        print("\nValidating with nucleotide pHMMS. Initiating nhmmscan.");log_file.write("\nValidating with nucleotide pHMMS. Initiating nhmmscan.")
        subprocess.call("nhmmscan -o /dev/null -E 0.00001 --cpu " + str(num_threads) + " --dfamtblout OUT_" + seq_file.replace(".fas", ".out") + " " + usr_args.dna_hmms + "nucl_pHMMs " + seq_file, shell=True)
        subprocess.call("grep -v \"^#\" OUT_" + (seq_file.replace(".fas", ".out")) + " > " + (seq_file.replace(".fas", "_nucleotide_pHMM.out")), shell=True)
        subprocess.call("rm OUT_" + (seq_file.replace(".fas", ".out")), shell=True)

        if (int(subprocess.check_output("ls -al " + (seq_file.replace(".fas", "_nucleotide_pHMM.out")) + " | awk \'{print $5}\'", shell=True)) > 0):
                total_ann = (str(subprocess.check_output("wc -l " + seq_file.replace(".fas", "_nucleotide_pHMM.out") + " | awk \'{print $1}\'", shell=True))).rstrip("\n")
                ltr_ann = (str(subprocess.check_output("egrep \"ERV|Long Terminal Repeat|LTR|HERV\" " + seq_file.replace(".fas", "_nucleotide_pHMM.out") + " | wc -l", shell=True))).rstrip("\n")
                if (ltr_ann == "0"):
                        best_hit = "NULL"
                        strand = "NULL"
                        e_val = "NULL"
                        hit_pos = "NULL"
                else:
                        top_ltr = (str(subprocess.check_output("egrep \"ERV|Long Terminal Repeat|LTR|HERV\" " + seq_file.replace(".fas", "_nucleotide_pHMM.out") + " | head -n 1", shell=True))).rstrip("\n")
                        best_hit = (str(subprocess.check_output("egrep \"ERV|Long Terminal Repeat|LTR|HERV\" " + seq_file.replace(".fas", "_nucleotide_pHMM.out") + " | head -n 1 | awk \'{print $1 \"-\" $15}\'", shell=True))).rstrip("\n")
                        strand = (str(subprocess.check_output("egrep \"ERV|Long Terminal Repeat|LTR|HERV\" " + seq_file.replace(".fas", "_nucleotide_pHMM.out") + " | head -n 1 | awk \'{print $5}\'", shell=True))).rstrip("\n")
                        e_val = (str(subprocess.check_output("egrep \"ERV|Long Terminal Repeat|LTR|HERV\" " + seq_file.replace(".fas", "_nucleotide_pHMM.out") + " | head -n 1 | awk \'{print $9}\'", shell=True))).rstrip("\n")
                        hit_pos = (str(((subprocess.check_output("grep -n \"" + top_ltr + "\" " + seq_file.replace(".fas", "_nucleotide_pHMM.out"), shell=True)).split(":"))[0])).rstrip("\n")
                nucl_file.write(total_ann + "\t" + ltr_ann + "\t" + best_hit + "\t" + str(e_val) + "\t" + strand + "\t" + str(hit_pos) + "\n")
        else:
                nucl_file.write("0\t0\tNULL\tNULL\tNULL\tNULL\n")

        
        ## tRNA HMM analysis
        print("Attempting to identify the PBS with tRNA pHMMS. Initiating nhmmscan.");log_file.write("Attempting to identify the PBS with tRNA pHMMS. Initiating nhmmscan.")
        subprocess.call("nhmmscan -o /dev/null -E 0.00001 --cpu " + str(num_threads) + " --dfamtblout OUT_" + seq_file.replace(".fas", ".out") + " " + usr_args.trna_hmms + "tRNA_pHMMs " + seq_file, shell=True)
        subprocess.call("grep -v \"^#\" OUT_" + (seq_file.replace(".fas", ".out")) + " > " + (seq_file.replace(".fas", "_tRNA_pHMM.out")), shell=True)
        subprocess.call("rm OUT_" + (seq_file.replace(".fas", ".out")), shell=True)

        if (int(subprocess.check_output("ls -al " + (seq_file.replace(".fas", "_tRNA_pHMM.out")) + " | awk \'{print $5}\'", shell=True)) > 0):
                res = ((subprocess.check_output("head -n 1 " + (seq_file.replace(".fas", "_tRNA_pHMM.out")) + " | awk \'{print $1}\'", shell=True)).rstrip("\n"))[:3]
                e_val = subprocess.check_output("head -n 1 " + (seq_file.replace(".fas", "_tRNA_pHMM.out")) + " | awk \'{print $5}\'", shell=True).rstrip("\n")
                strand = subprocess.check_output("head -n 1 " + (seq_file.replace(".fas", "_tRNA_pHMM.out")) + " | awk \'{print $6}\'", shell=True).rstrip("\n")
                tRNA_file.write(str(res) + "\t" + str(e_val) + "\t" + str(strand) + "\n")
        else:
                tRNA_file.write("NULL\tNULL\tNULL\n")
                
        ## Protein HMM analysis
        # calls the EMBOSS tool 'transeq' to translate the sequence into all 6 aa frames. The -clean flag changes all stops (*) to 'unknown' (X) as some programs don't like the presence of stop codons
        print("Validating with protein pHMMs. Translating " + seq_file + " using EMBOSS transeq protocol.");log_file.write("Validating with protein pHMMs. Translating " + seq_file + " using EMBOSS transeq protocol.\n")
        subprocess.call("transeq -sequence " + seq_file + " -frame 6 -clean Y -outseq translated_seq.fas", shell=True)
        print("Initiating hmmscan to scan translated sequences for protein profiles defined in the \"prot_pHMMs\" database.");log_file.write("Initiating hmmscan to scan translated sequences for protein profiles defined in the \"prot_HMMs\" database.")
        subprocess.call("hmmscan -o /dev/null -E 0.00001 --cpu " + str(num_threads) + " --domtblout OUT_" + seq_file.replace(".fas", ".out") + " " + usr_args.protein_hmms + "prot_pHMMs translated_seq.fas", shell=True)
        subprocess.call("grep -v \"^#\" OUT_" + (seq_file.replace(".fas", ".out")) + " | sort -g -k7 > " + (seq_file.replace(".fas", "_protein_pHMM.out")), shell=True)
        subprocess.call("rm OUT_" + (seq_file.replace(".fas", ".out")), shell=True)

        domain_groups = ["GAG", "AP", "INT|GIN", "RT", "RNaseH", "DUT", "ENV", "ATF", "CHR", "MOV", "TAV", "VAP", "BEL|NEF|ORF|REV|REX|ROF|SORF|TAT|TAX|TOF|VIF|VPR|VPX"]
        if (int(subprocess.check_output("ls -al " + (seq_file.replace(".fas", "_protein_pHMM.out")) + " | awk \'{print $5}\'", shell=True)) > 0):
                full_report = ""
                i = 1
                for domain in domain_groups:
                        if (i < len(domain_groups)):
                                full_report += prot_domain_sorter(domain, (seq_file.replace(".fas", "_protein_pHMM.out")), False)
                        else:
                                full_report += prot_domain_sorter(domain, (seq_file.replace(".fas", "_protein_pHMM.out")), True)
                        i+=1
                prot_file.write(full_report.rstrip("\t") + "\n")             
        else:
                prot_file.write("NULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n")

        ## RepeatMasker analysis
        # calls RepeatMasker on each sequence and extracts the annotations with total length covered
        print("Analysing sequence against RepBase motif databases.");log_file.write("Analysing sequence against RepBase motif databases.")
        subprocess.call("repeatmasker -species " + usr_args.species + " " + seq_file + " &> /dev/null", shell=True)
        
        if os.path.isfile("./" + seq_file + ".out"):
                subprocess.call("tail -n +4 " + seq_file + ".out | awk \'{print $11 \"\t\" (($7+1)-$6) \"\t\" $9}\' | sort > RM-OUT.out", shell=True)
                subprocess.call("rm " + seq_file + ".*", shell=True)
                all_elements = []
                for line in (open("RM-OUT.out").read().rstrip("\n").split("\n")):
                        elements = line.split("\t")
                        all_elements += elements

                z=[]
                y = []
                i=0
                j=0
                k=0
                rm_strand = ""
                for x in all_elements:
                        x = x.replace("?", "")
                        if x.isdigit():
                                j += int(x)
                                k+=1
                        elif (x=='+' or x=='C'):
                                rm_strand += str(x)
                        else:
                                rep_class = x
                                if re.search(r"/", x):
                                        rep_class = (x.split("/"))[0]
		
                                if (i==0):
                                        if re.match("DNA", rep_class):
                                        	z += ["DNA_transposon"]
                                        else:
                                                z += [rep_class]
		
                                if ((i>= 3) and (rep_class != ((all_elements[i-3]).split("/")[0]))):
                                        cov = str(round(((100 / seq_length) * int(j)), 2)) + "%"
                                        z += [int(j)]
                                        z += [cov]
                                        z += [int(k)]
                                        z += [rm_strand]
                                        y += [z]
                                        z=[]
                                        j=0
                                        k=0
                                        if re.match("DNA", rep_class):
                                        	z += ["DNA_transposon"]
                                        else:
                                                z += [rep_class]
                        if (i==(len(all_elements) -1)):
                                cov = str(round(((100 / seq_length) * int(j)), 2)) + "%"
                                z += [int(j)]
                                z += [cov]
                                z += [int(k)]
                                z += [rm_strand]
                                y += [z]		
                        i+=1	
                y.sort(lambda l1, l2: cmp(l1[1], l2[1]) or (l1[1] == l2[1] and cmp(l1[0], l2[0])))
                processed_list = y[::-1]
                rm_results = ""
                ltr_percent = "0"
                ltr_num = "NULL"
                ltr_strand = "NULL"
                i=1
                for element in processed_list:
                        if (element[0] == "LTR"):
                                ltr_percent = (element[2]).replace("%", "")
                                ltr_num = str(i)
                                #get ltr strand
                                plus_count = (element[4]).count("+")
                                C_count = (element[4]).count("C")
                                strand_length = plus_count + C_count
                                if plus_count >= C_count:
                                        ltr_strand = "+"
                                else:
                                        ltr_strand = "-"
                        rm_results += str(element).replace("'", "").replace("[", "").replace("]", "")
                        if (i < len(processed_list)):
                                rm_results += "; "
                        i+=1

                rm_file.write(rm_results + "\t" + ltr_percent + "\t" + ltr_num + "\t" + ltr_strand + "\n")
                subprocess.call("rm RM-OUT.out", shell=True)
        else:
                rm_file.write("NULL\tNULL\tNULL\tNULL\n")

print("\nValidation complete. Formatting results file."); log_file.write("\nValidation complete. Formatting results file.")

for file_name in os.listdir("."):
        if (file_name.endswith(".fas") or file_name.endswith("_pHMM.out")):
                os.remove(file_name)
if os.path.isfile("5_col_pos_file.txt"):
	os.remove("5_col_pos_file.txt")

seq_info_file.close()
rm_file.close()
nucl_file.close()
tRNA_file.close()
prot_file.close()
subprocess.call("paste seq_info.txt repeatmasker_scores.txt validated_nucl_scores.txt validated_tRNA_scores.txt validated_prot_scores.txt > domain_validated_scores.txt", shell=True)

to_remove = ""
with open("domain_validated_scores.txt") as f:
        reader = csv.reader(f, delimiter="\t")
        cols_to_remove = []
        i=1
        for col in zip(*reader):
                column = list(OrderedDict.fromkeys((list(col)[2:])))
                column.sort
		if i>15:                
			if ((len(column) == 2) and (((column[0]== "0") and (column[1] == "NULL")) or ((column[1]== "0") and (column[0] == "NULL")))):
				cols_to_remove += [i]
			if (len(column) == 1) and (column[0]== "0" or column[0] == "NULL"):
				cols_to_remove += [i]
		i+=1
	to_remove = (str(cols_to_remove)).replace("[", "").replace("]", "").replace(" ", "")

subprocess.call("cut --complement -f" + to_remove + " domain_validated_scores.txt > domain_scores.txt", shell=True)

print("\nPerforming automated analysis of validation output. Result in column 1 of \"domain_validated_hits.txt\".");log_file.write("\nPerforming automated analysis of validation output. Result in column 1 of \"domain_validated_hits.txt\".")
valid_file = open("valid_scores.txt", "w")
valid_file.write("-----\nValidated?\n")
validation_file = (open("domain_scores.txt").read().rstrip("\n")).split("\n")

# check if tRNA columns have been removed from file due to having no matches
tRNA=True
if (((validation_file[0]).split("\t"))[16] != "--tRNA--"):
        tRNA=False
validation_file.pop(0)
validation_file.pop(0)

for line in validation_file:
        line_split = line.split("\t")
        R = False
        D = False
        P = False
        t = False
        D_pos = ""

        if (line_split[8] == "1"):
                R=True

        if ((line_split[11]) != "NULL"):
		if (int(line_split[11]) > 0):
                	D = True
                	D_pos = str(line_split[15])
	
        if (len(line_split) > 17):
                if (tRNA==True):
                        if (line_split[16] != "NULL"):
                                t = True
                        if (line_split[19] != "NULL"):
                                P = True
                else:
                        if (line_split[16] != "NULL"):
                                P = True

        if (R==True and D==True and P==True and t==True):
                valid_file.write("D" + D_pos + "+P+R+t\n")
        elif (R==True and D==True and P==True and t==False):
                valid_file.write("D" + D_pos + "+P+R\n")
        elif (R==True and D==True and P==False and t==True):
                valid_file.write("D" + D_pos + "+R+t\n")
        elif (R==True and D==False and P==True and t==True):
                valid_file.write("P+R+t\n")
        elif (R==True and D==True and P==False and t==False):
                valid_file.write("D" + D_pos + "+R\n")
        elif (R==True and D==False and P==True and t==False):
                valid_file.write("P+R\n")
        elif (R==True and D==False and P==False and t==False):
                valid_file.write("R\n")
        elif (R==True and D==False and P==False and t==True):
                valid_file.write("R+t\n")
        elif (R==False and D==True and P==True and t==True):
                valid_file.write("D" + D_pos + "+P+t\n")
        elif (R==False and D==True and P==True and t==False):
                valid_file.write("D" + D_pos + "+P\n")
        elif (R==False and D==True and P==False and t==True):
                valid_file.write("D" + D_pos + "+t\n")
        elif (R==False and D==False and P==True and t==True):
                valid_file.write("P+t\n")
        elif (R==False and D==True and P==False and t==False):
                valid_file.write("D" + D_pos + "\n")
        elif (R==False and D==False and P==True and t==False):
                valid_file.write("P\n")
        elif (R==False and D==False and P==False and t==True):
                valid_file.write("t\n")
        else:
                valid_file.write("XX\n")                
valid_file.close()

subprocess.call("paste valid_scores.txt domain_scores.txt > " + usr_args.out_prefix + "_full_output.txt", shell=True)
subprocess.call("grep -E -v \"^XX|^t\" " + usr_args.out_prefix + "_full_output.txt > " + usr_args.out_prefix + "_potentially_validated_output.txt", shell=True)
subprocess.call("rm *_scores.txt seq_info.txt", shell=True)

end = subprocess.check_output("date | awk \'{print $3 \"th \" $2 \" \" $6 \"  \"  $4}\'", shell=True)
print("\nProcess complete.\nProcess began: " + str(startdate) + "Process finished: " + str(end) + "\n");log_file.write("\nProcess complete.\nProcess began: " + str(startdate) + "Process finished: " + str(end) + "\n")
log_file.close()