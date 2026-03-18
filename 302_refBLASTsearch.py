### python 302_refBLASTsearch.py [-h] [--threads THREADS] /full/path/reference_multi_seq.fasta /full/path/reference_genome.fasta &> refBLAST.log &

## defining user argument parsing

from __future__ import division
import argparse
parser = argparse.ArgumentParser(
	description="refBLASTsearch takes the reference sequences and searches the given genome for all homologous hits.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ref_seq", help="Specify the reference multi sequence fastA file - use full path, not relative")
parser.add_argument("ref_genome", help="Specify your reference organism genome file - use full path, not relative")
parser.add_argument("--threads", type=int, default=8, help="Specify the number of threads (default = 8)")

skip_choice = parser.add_mutually_exclusive_group()
skip_choice.add_argument("--skip_ref_blast", action="store_true", help="skip the initial BLAST - choose one skip only")
skip_choice.add_argument("--skip_pos_merge", action="store_true", help="skip the initial BLAST AND the pos_merge protocol - choose one skip only")
skip_choice.add_argument("--skip_validate_by_homology", action="store_true", help="skip the initial BLAST AND the pos_merge protocol AND re-BLAST validation - choose one skip only")
skip_choice.add_argument("--skip_reduce_to_best", action="store_true", help="skip the initial BLAST AND the pos_merge protocol AND re-BLAST validation AND reduction to best hits/none matched - choose one skip only")
skip_choice.add_argument("--skip_best_rm", action="store_true", help="skip the initial BLAST AND the pos_merge protocol AND re-BLAST validation AND reduction to best hits/none matched AND best/none RM - choose one skip only")

usr_args = parser.parse_args()

import os
import re
import subprocess
import sys
import time

os.mkdir("refBLAST_files")
os.chdir("./refBLAST_files")

startdate = subprocess.check_output("date | awk \'{print $3 \" \" $2 \" \" $6 \"  \"  $4}\'", shell=True)

sys.path.append("full_path_to_LocaTR/static_functions.py")
import static_functions

######### Pre-process checks and initiation of "log_file.txt" #########

if os.path.isfile("log_file.txt"):                      # log_file.txt initiation
	log_file = open("log_file.txt", "a")
	log_file.write("\n\n\n\n\n------------------------------------------------------------------------------------------\n\n\n" + str(startdate))
else:
	log_file = open("log_file.txt", "w")
	log_file.write(str(startdate))

# Check if reference genome exists as BLAST database already. If not, create BLAST database files.
if os.path.isfile(usr_args.ref_genome + ".nhr") and os.path.isfile(usr_args.ref_genome + ".nin") and os.path.isfile(usr_args.ref_genome + ".nsq"):
	print("\nAll BLAST database files found. Continuing to create query files.");log_file.write("\nAll BLAST database files found. Continuing to create query files.")	
else:
	print("\nBLAST database files not present. Implementing makeblastdb protocol.");log_file.write("\nBLAST database files not present. Implementing makeblastdb protocol.")
	subprocess.call(("makeblastdb -in " + usr_args.ref_genome + " -dbtype nucl"), shell=True)
	print("BLAST database files created. Continuing to create query files.");log_file.write("BLAST database files created. Continuing to create query files.")


######### Initial Definitions used throughout process #########
print("\nFunctions defined. Defining key variables.");log_file.write("\nFunctions defined. Defining key variables.")

ref_seq_file = open(usr_args.ref_seq, "rU").read().rstrip("\n")      # define and read the reference sequences file - used as queries in the initial search and references in validation
ref_gen_file = open(usr_args.ref_genome).read().rstrip("\n")            # define and read the reference genome file

num_processes = usr_args.threads                   # defines the number of processes available to BLAST during the initial search and validation

# Create reference genome dict for sequence extraction
print("Indexing reference genome for rapid sequence extraction.");log_file.write("Indexing reference genome for rapid sequence extraction.")
ref_gen_headers = static_functions.header_extractor(ref_gen_file)
ref_gen_seq_list = static_functions.seq_only_extractor(ref_gen_file)
ref_gen_dict = static_functions.seq_dict_creator(usr_args.ref_genome)

### SKIP 1
if (usr_args.skip_ref_blast or usr_args.skip_pos_merge or usr_args.skip_validate_by_homology or usr_args.skip_reduce_to_best or usr_args.skip_best_rm):
        print("\nUser has stated to skip the initial BLAST protocol. Moving on to pos_merge protocol.");log_file.write("\nUser has stated to skip the initial BLAST protocol. Moving on to pos_merge protocol.\n")
else:

        # Create separate files for each of the reference sequence fasta sequences
        print("Indexing of reference genome complete. Creating individual query files for each reference sequence.");log_file.write("Indexing of reference genome complete. Creating individual query files for each reference sequence.")
        ref_seq_files = static_functions.fasta_splitter(ref_seq_file, "Y")
        print(str(len(ref_seq_files)) + " query files ready for BLAST analysis.\n");log_file.write(str(len(ref_seq_files)) + " query files ready for BLAST analysis.\n")

        ######### refBLASTsearch protocol #########

        ## PART 1 - Perform BLAST and create positions file

        print("Initiating refBLASTsearch protocol with blastn and tblastx algorithms.\n");log_file.write("Initiating refBLASTsearch protocol with blastn and tblastx algorithms.\n")
        methods = ["blastn", "tblastx"]
        pos_file = open("refBLAST_positions.txt", "w")

        i=1
        for query in ref_seq_files:
                s = time.time()
                print("Loading query file " + str(i) + " of " + str(len(ref_seq_files)) + ".");log_file.write("Loading query file " + str(i) + " of " + str(len(ref_seq_files)) + ". ")
                i+=1
                # perform both blastn and tblastx for each query file - use 10^-10 E value threshold and store contig, start, end and length of match
                for method in methods:
                        print("Initiating " + method + "...");log_file.write("Initiating " + method + "...")
                        subprocess.call(method + " -db " + (usr_args.ref_genome) + " -query " + query + " -evalue 0.0000000001 -num_threads " + str(num_processes) + " -out results.out -outfmt \"6 sacc sstart send length\" -max_target_seqs 1000000", shell=True)
                        hits = (open("results.out").read().rstrip("\n").split("\n"))
                        # if any matches are returned for each hit remove any matches shorter than 100bp and apply strandedness based on whether start is smaller than end (+) or not (-). Write matches to file
                        if hits[0] != "":
                                for hit in hits:
                                        cols = hit.split("\t")
                                        if cols[3] >= 100:
                                                if (cols[1] < cols[2]):
                                                        pos_file.write(str(cols[0]) + "\t" + str(cols[1]) + "\t" + str(cols[2]) + "\t+\tBR\n")
                                                else:
                                                        pos_file.write(str(cols[0]) + "\t" + str(cols[2]) + "\t" + str(cols[1]) + "\t-\tBR\n")
                        print("..." + method + " search completed.");log_file.write("..." + method + " search completed. ")
                print("Query processing took " + "{0:.3f}".format(time.time() - s) + " seconds.\n");log_file.write("Query processing took " + "{0:.3f}".format(time.time() - s) + " seconds.\n")
	
        pos_file.close()
        os.remove("results.out")

### END OF SKIP 1
### SKIP 2 
if (usr_args.skip_pos_merge or usr_args.skip_validate_by_homology or usr_args.skip_reduce_to_best or usr_args.skip_best_rm):
        print("User has stated to skip the pos_merge protocol. Moving on to validate_by_homology.");log_file.write("\nUser has stated to skip the pos_merge protocol. Moving on to validate_by_homology.\n")
else:


        ## PART 2 - Merge positions found by BLAST searches. Merging protocol divided by chromosome for more efficient sorting.

        # If a RM analysis has been completed to supplement the reference sequences it gets added in here. "rm_hits.txt" is the output of "run_repeatmasker.py". Else, the output from above is renamed to match.
        if os.path.isfile("./rm_hits.txt"):
                subprocess.call("cat refBLAST_positions.txt rm_hits.txt > RB_plus_RM.txt", shell=True)
        else:
                subprocess.call("cp refBLAST_positions.txt RB_plus_RM.txt", shell=True)
        #format the list using the static_functions defined list_formatter() function once the file has been sorted and duplicated lines removed - list is split by contig/chromosome
        pos_by_chro = static_functions.list_formatter((subprocess.check_output("sort -k1 -n RB_plus_RM.txt | uniq", shell=True)), 0)
        length = (str(subprocess.check_output("wc -l refBLAST_positions.txt | awk \'{print $1}\'", shell=True))).rstrip("\n")
        print("\nBLAST searches complete.\nThere are " + length + " results in the positions file.\nOriginal positions list has been sorted ready for merging. Proceeding with positions_merger protocol.\n");log_file.write("BLAST searches complete.\nThere are " + length + " results in the positions file.\nOriginal positions list has been sorted ready for merging. Proceeding with positions_merger protocol.\n")

        all_positions = []
        # for each contig/chromosome with results, if there are more than 1000 hits, split into sublists of max 1000 hits and perform pos_merge protocol per sublist, then merge sublists with pos_merge to sort boundaries
        # if there are 1000 or fewer hits do one pos_merge protocol. All pos then concatenated together and written to file.
        # pos_merge i)removes positions completely contained in other positions, or ii) merges positions which overlap, touch or are separated by less than 11bp. The strandedness is resolved by going with the most common, if equal + is taken as default
        for x in pos_by_chro:
                if len(x) > 1000:
                        interim_pos = []
                        x_split = [x[j:j+1000] for j in range(0, len(x), 1000)]
                        print("Processing contig \"" + str(x[0][0]) + "\" - " + str(len(x)) + " positions to process, split across " + str(len(x_split)) + " sublists.")
                        log_file.write("Processing contig \"" + str(x[0][0]) + "\" - " + str(len(x)) + " positions to process, split across " + str(len(x_split)) + " sublists.")
                        m=1
                        for i in x_split:
                                interim_pos += static_functions.positions_merger(i, log_file)
                                print("  Processed sublist " + str(m) + " of " + str(len(x_split)))
                                log_file.write("  Processed sublist " + str(m) + " of " + str(len(x_split)))
                                m+=1
                        print("  \"" + str(x[0][0]) + "\" sublists concatenated. Merging sublists.\n  Original list length: " + str(len(x)) + "\n  Current list length: " + str(len(interim_pos)))
                        log_file.write("  \"" + str(x[0][0]) + "\" sublists concatenated. Merging sublists.\n  Original list length: " + str(len(x)) + "\n  Current list length: " + str(len(interim_pos)))
                        all_pos = static_functions.positions_merger(interim_pos, log_file)
                        all_positions += all_pos
                        print("  Final list length: " + str(len(all_pos)));log_file.write("  Final list length: " + str(len(all_pos)))
                        print("  Contig merge complete.\n");log_file.write("  Contig merge complete.\n")
                else:
                        print("Processing contig \"" + str(x[0][0]) + "\" - " + str(len(x)) + " positions to process.")
                        log_file.write("Processing contig \"" + str(x[0][0]) + "\" - " + str(len(x)) + " positions to process.")
                        all_pos = static_functions.positions_merger(x, log_file)
                        all_positions += all_pos
                        print("  Final list length: " + str(len(all_pos)));log_file.write("  Final list length: " + str(len(all_pos)))
                        print("  Contig merge complete.\n");log_file.write("  Contig merge complete.\n")
        print("Merge protocol complete. Writing positions to file.\n");log_file.write("Merge protocol complete. Writing positions to file.\n")

        outfile = open("refBLAST_merged_positions.txt", "w")
        for pos in all_positions:
                outfile.write(str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + str(pos[3]) + "\t" + str(pos[4]) +"\n")
        outfile.close()


        ## PART 3 - Extract refBLAST_merged_positions sequences from the reference genome

        print("Extracting putative sequences from the reference genome.");log_file.write("Extracting putative sequences from the reference genome. ")
        print("Genomic sequence loaded. Position information formatted. Writing putative sequences to file.");log_file.write("Genomic sequence loaded. Position information formatted. Writing putative sequences to file.")
        static_functions.seq_extract("refBLAST_seq.fasta", all_positions, ref_gen_dict, "Y")
        print("Sequences written to file.\n");log_file.write("Sequences written to file.\n")


### END OF SKIP 2
### SKIP 3
if (usr_args.skip_validate_by_homology or usr_args.skip_reduce_to_best or usr_args.skip_best_rm):
        print("User has stated to skip the validate_to_homology protocol. Moving on to reduce_to_best protocol.");log_file.write("\nUser has stated to skip the validate_to_homology protocol. Moving on to reduce_to_best protocol.\n")
else:

        ######### validate_by_homology protocol #########

        ## PART 1 - Perform tBLASTx on all reference query files against a refBLAST_seq BLAST database

        # create a BLAST database from identified, merged sequences. use this for reciprocal BLAST using the original reference sequences - are all putative positions retrieved?
        print("\nInitiating validate_by_homology protocol on extracted sequences.\n");log_file.write("Initiating validate_by_homology protocol on extracted sequences.\n")
        print("Creating BLAST database from putative LTR retrotransposon sequences.");log_file.write("Creating BLAST database from putative LTR retrotransposon sequences.")
        subprocess.call(("makeblastdb -in refBLAST_seq.fasta -dbtype nucl"), shell=True)
        print("\nBLAST database created.");log_file.write("\nBLAST database created.")

        hit_names = static_functions.header_extractor(open("refBLAST_seq.fasta").read().rstrip("\n"))
        print("Query files created. There are " + str(len(ref_seq_files)) + " query sequences.\nImplementing tBLASTx using all query files against the putative element sequence BLAST database.\n");log_file.write("Query files created. There are " + str(len(ref_seq_files)) + " query sequences.\nImplementing tBLASTx using all query files against the putative element sequence BLAST database.\n")

        j=0
        out_names = []
        for filename in ref_seq_files:
                core_name = re.sub(".fas", "", filename)
                results_file = "results_" + core_name + ".out"
                # perform tblastx only using E value 10^-10, 8 threads, allowing for 5000000 matches, creating output file of query name etc.
                subprocess.call(("tblastx -db refBLAST_seq.fasta -query " + filename + " -evalue 0.0000000001 -num_threads 8 -out " + results_file + " -outfmt \"6 qacc sacc qstart qend sstart send evalue length\" -max_target_seqs 5000000"), shell=True)
                # remove hits shorter than 30bp and then create list with subject name, e value, start and end
                filt_hits = ((subprocess.check_output(("sort -t 2 -n " + results_file + " | awk \'$8>30{print $2 \"\t\" $7 \"\t\" $5 \"\t\" $6}\'"), shell=True)).rstrip("\n")).split("\n")
                subprocess.call(("rm " + filename + " " + results_file), shell=True)
	
                if filt_hits[0] != '':
                        blast_hit_list = []
                        for hit in filt_hits:
                                blast_hit_list += [(hit.split("\t"))]

                        # for each hit create a dict with its E value, start and end. If the evalue is 0.0 (aka exact match), this is converted to 1e-200 for later use
                        hit_dict = {}
                        for hit in blast_hit_list:
                                homology = ""
                                if (str(hit[1]) == "0.0"):
                                        homology = "1e-200"
                                else:
                                        homology = str(hit[1])
                                hit_dict[(hit[0])] = str(homology + "\t" + hit[2] + "\t" + hit[3])
	
                # each identified position is then scored against each ref sequence in this manner, creating a full list of scores. If there is no match (none at all or for one sequence) you give a blank score with an E value of 100
                out = "formatted_" + core_name + ".out"
                out_names += [out]
                final_output = open((out), "w")
                for hit_name in hit_names:
                        if (filt_hits[0] == '') or (hit_dict.get(hit_name) == None):
                                final_output.write(str(100) + "\t0\t0\n")	
                        else:
                                final_output.write(str(hit_dict[hit_name]) + "\n")	
                final_output.close()
                j+=1
	
                if (j%10) == 0:
                        print("Processed " + str(j) + " query files of " + str(len(ref_seq_files)) + ".\n");log_file.write("Processed " + str(j) + " query files of " + str(len(ref_seq_files)) + ".\n")

        # create file with all sequence names (identified in the refBLAST)
        print("Writing results to file.");log_file.write("Writing results to file.")
        first_col = open("names.out", "w")
        for name in hit_names:
        	first_col.write(name + "\n")
        first_col.close()

        # create a top line of appropriate headers. Name-E is the evalue, -sub_s is subject start pos, -sub_e is subject end pos
        top_line = "Sequence"
        for names in ref_seq_files:
                name = re.sub(".fas", "", names)
        	top_line += "\t" + name + "-E\t" + name + "-sub_s\t" + name + "-sub_e"
        top = open("top.out", "w")
        top.write(top_line + "\n")
        top.close()

        # paste all files giving results for each reference sequence
        i=0
        for filename in out_names:
                if i==0:
                        subprocess.call("paste names.out " + filename + " > processed.tmp", shell=True)
                else:
                        subprocess.call("paste processed.tmp " + filename + " > processed.out", shell=True)
                        subprocess.call("mv processed.out processed.tmp", shell=True)
                i+=1

        # add in the top line of headers and do some file removal
        subprocess.call("cat top.out processed.tmp > processed_scores.txt", shell=True)
        subprocess.call("rm *.out *.fasta.n* *.tmp", shell=True)
        print("Full homology validation scores produced and written to file.");log_file.write(" Full homology validation scores produced and written to file. ")

### END OF SKIP 3
### SKIP 4
if (usr_args.skip_reduce_to_best or usr_args.skip_best_rm):
        print("User has stated to skip the reduce_to_best protocol. Moving on to element validation.");log_file.write("\nUser has stated to skip the reduce_to_best protocol. Moving on to element validation.\n")
else:

        ## PART 2 - Extract information on best hits for each fully validated putative element

        print("\nReducing results to best matches only.\n");log_file.write("Reducing results to best matches only.\n")
        results = (open("processed_scores.txt").read().rstrip("\n")).split("\n")
        Top = (results[0]).split("\t")
        Top.pop(0)
        # extract all headers from the processed_scores.txt file, remove the header name 'Sequences' and then extract just the ref header names (1 of each, not the three copies)
        topline = []
        topline.extend(Top[0::3])
        ref_names = []
        for element in topline:
        	ref_names += [element.replace(".out-E", "")]
        # remove the first column of sequence names
        results.pop(0)

        no_match = open("none_matched.txt", "w")
        seq_with_matches = []
        for row in results:
                line = row.split("\t")
                # applicable replace for one of the linkage groups in galGal4
                seq = (line[0]).replace("28_E5", "28.E5")
                line.pop(0)
                each_query = [line[i:i+3] for i in range(0, len(line), 3)]
                i=0
                matches = []
                for query in each_query:
                        if query[0] != '100':
                                query.insert(0, ref_names[i])
                                matches += [(query)]
                        i+=1
	
                contig = ((seq.split("_"))[0]).replace(".", "_")
                startpos = (((seq.split("_"))[1]).split(":"))[0]
                endpos = (((seq.split("_"))[1]).split(":"))[1]
                strand = (seq.split("_"))[2]
	
                seq_header = str(contig) + "\t" + str(startpos) + "\t" + str(endpos) + "\t" + strand
	
                if matches:
                        matches.sort(lambda l1, l2: cmp(float(l1[1]), float(l2[1]) or (float(l1[1]) == float(l2[1]) and cmp(l1[0], l2[0]))))
                        matches.insert(0, seq_header)
                        seq_length = abs(int((((seq.split(":"))[1]).split("_"))[0]) - int((((seq.split(":"))[0]).split("_"))[-1])) + 1
                        matches.insert(1, seq_length)
		
                        coverage = 0
                        seq_range = ""
                        best_evalue = ""
                        evalue = []
		
                        if len(matches) < 4:
                                coverage = round(((100 / seq_length) * ((abs(int(matches[2][3]) - int(matches[2][2])) + 1))), 2)
			
                                if int(matches[2][3]) < int(matches[2][2]):
                                        st = int(matches[2][3])
                                        end = int(matches[2][2])
                                else:
                                        st = int(matches[2][2])
                                        end = int(matches[2][3])
                                seq_range = str(st) + "-" + str(end)
                                
                                best_evalue = str((matches[2][0]).replace("_", "-")) + "_" + str(matches[2][1]) + "_" + str(matches[2][2]) + ":" + str(matches[2][3])
                                matches.insert(2, coverage)
                                matches.insert(3, seq_range)
                                matches.insert(4, best_evalue)
		
                        else:
                                cov = []
                                matches_to_sort = matches[2:]
                                best_evalue = str((matches_to_sort[0][0]).replace("_", "-")) + "_" + str(matches_to_sort[0][1]) + "_" + str(matches_to_sort[0][2]) + ":" + str(matches_to_sort[0][3])
                                x=0
                                for pos in matches_to_sort:
                                        if int(pos[2]) > int(pos[3]):
                                                st, end = 2, 3
                                                pos[end], pos[st] = pos[st], pos[end]
                                        cov += [[int(pos[2]),int(pos[3]),int(x)]]
                                        x+=1
			
                                cov.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]))))
                                sort=None
                                z=0
                                while z < (len(cov) - 1):
                                        if sort == True:
                                                cov.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]))))
				
                                        if (int(cov[z+1][1]) <= int(cov[z][1])):
                                                cov.pop(z+1)
                                                sort=True
					
                                        elif (int(cov[z][1]) >= int(cov[z+1][0])) or (int(cov[z+1][0]) == (int(cov[z][1] + 1))):
                                                cov += [[int(cov[z][0]),int(cov[z+1][1])]]
                                                cov.pop(z+1)
                                                cov.pop(z)
                                                sort=True
                                        else:
                                                z+=1
                                                sort=False
                                total_cov = 0
                                seq_coverage = ""
                                for pos in cov:
                                        total_cov += (int(pos[1]) - int(pos[0])) + 1
                                        seq_coverage += str(pos[0]) + "-" + str(pos[1]) + ","
                                coverage = round(((100 / seq_length) * total_cov), 2)
                                seq_range = re.sub(",$", "", seq_coverage)
			
                                matches.insert(2, coverage)
                                matches.insert(3, seq_range)
                                matches.insert(4, best_evalue)

                        seq_with_matches += [matches]
		
                else:
                        no_match.write(str(seq_header) + "\trB_hom-NM\n")

        no_match.close()
	
        k=0
        y=0
        for m in seq_with_matches:
                if (m[4]).startswith("CR1"):
                        seq_with_matches.pop(k)
                        y+=1
                k+=1
	
        print(str(len(hit_names)) + " sequences were used to build the BLAST database.\n" + str(len(seq_with_matches)) + " have been matched by reference sequences in this validation.\n" + str(y) + " matches have been removed due to primary homology with CR1 elements.\n\nInformation for " + str(len(seq_with_matches) - y) + " sequences being written to file.\n");log_file.write(str(len(hit_names)) + " sequences were used to build the BLAST database.\n" + str(len(seq_with_matches)) + " have been matched by reference sequences in this validation.\n" + str(y) + " matches have been removed due to primary homology with CR1 elements.\n\nInformation for " + str(len(seq_with_matches) - y) + " sequences being written to file.\n")	

        output = open("best_hits.txt", "w")
        for line in seq_with_matches:
                line_len = len(line)
                i=0
                line_to_write = ""
                for element in line:
                        if i == len(line) - 1:
                                line_to_write += (str(line[i]) + "\n")
                        else:
                                line_to_write += (str(line[i]) + "\t")			
                        i+=1
                output.write(str(line_to_write))
        output.close()

### END OF SKIP 4
### SKIP 5
if (usr_args.skip_best_rm):
        print("User has stated to skip the best/none RM. Moving on to element validation.");log_file.write("\nUser has stated to skip the best/none RM. Moving on to element validation.\n")
else:

        ## PART 3 - extract sequences in files that were matched and not matched

        print("\nExtracting sequences for positions in \"best_hits.txt\" and \"none_matched.txt\".\n");log_file.write("Extracting sequences for positions in \"best_hits.txt\" and \"none_matched.txt\".\n")

        files = ["best_hits.txt", "none_matched.txt"]
        for start_file in files:
                subprocess.call("sed s/_1/.1/g " + start_file + " > sed_replaced_file.txt", shell=True)
                file_core = (start_file.split("."))[0]
                print("Extracting sequences for \"" + start_file + "\" to \"" + file_core + "_seq.fasta\".");log_file.write("Extracting sequences for \"" + start_file + "\" to \"" + file_core + "_seq.fasta\".")
                #opens existing positions file and processes into a list of positions
                pos_file = open("sed_replaced_file.txt").read().rstrip("\n")
                positions_split = (pos_file.replace("\t", ", ")).split("\n")
                pos_list = []
                for element in positions_split:
                        x = (element.split(", "))
                        positions = []
                        count = 0
                        for val in x:
                                value = val.replace("\r", "")
                                count = count + 1
                                if count == 2 or count == 3:
                                        y = int(value)
                                        positions = positions + [y]
                                #not changing to int for pos[0] or pos[3] which may use non number characters
                                else:
                                        positions = positions + [value]
                        pos_list = pos_list + [positions]
                static_functions.seq_extract(("./" + file_core + "_seq.fasta"), pos_list, ref_gen_dict, "Y")
                subprocess.call("rm sed_replaced_file.txt", shell=True)
                # run repeatmasker on file
                print("Running RepeatMasker on \"" + file_core + "_seq.fasta\".");log_file.write(" Running RepeatMasker on \"" + file_core + "_seq.fasta\". ")
                subprocess.call(("repeatmasker -species vertebrates -nolow " + file_core + "_seq.fasta"), shell=True)
                rm_files = [(file_core + "_seq.fasta.cat"),(file_core + "_seq.fasta.cat.gz"),(file_core + "_seq.fasta.masked"), (file_core + "_custom_seq.fasta.out"), (file_core + "_def_seq.fasta.out")]
                for rm_file in rm_files:
                        if os.path.isfile(rm_file):
                                subprocess.call("rm " + rm_file, shell=True)
                print("RepeatMasker complete.\n");log_file.write("RepeatMasker complete.\n")
        print("Sequence extract complete.\n");log_file.write("Sequence extract complete.\n")

### END OF SKIP 5

## PART 4 - Processing RepeatMasker output from the best matches and no matches runs to remove CR1 annotations and (generally short) reincoporate LTR annotations

# pre-processing for best_hits_seq.fasta
print("Processing \"best_hits_seq.fasta.out\" before removing CR1 annotated fragments.");log_file.write("Processing \"best_hits_seq.fasta.out\" before removing CR1 annotated fragments.") 
if os.path.isfile("CR1_hits.tmp"):
	subprocess.call("rm CR1_hits.tmp", shell=True)

subprocess.call("awk \'{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\tBR\"}\' best_hits.txt | sed s/_1/.1/g > best_hits_no_detail.txt", shell=True)
if int((subprocess.check_output("grep -m 1 \"CR1\" best_hits_seq.fasta.out | wc -l", shell=True)).rstrip("\n")) > 0:
        # grep all lines containing CR1 annotation and acquire a list of the containing fragment. Then grep each line from that sequence fragment and append to file.
        subprocess.call("for i in `grep \"CR1\" best_hits_seq.fasta.out | awk \'{print $5}\' | sort -n | uniq -c | awk \'{print $2}\'`;do grep $i best_hits_seq.fasta.out | awk \'{print $1\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$10\"\\t\"$11}\'>> CR1_hits.tmp;echo \"\" >> CR1_hits.tmp;done", shell=True)
        subprocess.call("cp CR1_hits.tmp best_seq_CR1_annotations.txt", shell=True)
        static_functions.cr1_fragment_remover("best_seq_CR1_annotations.txt", "best_hits_no_detail.txt", "01", log_file, "BR")
else:
        print("No CR1 annotations found. Moving on to \"none_matched_seq.fasta.out\" analysis.\n");log_file.write("No CR1 annotations found. Moving on to \"none_matched_seq.fasta.out\" analysis.\n")

#pre-processing for none_matched.fasta
print("Processing \"none_matched_seq.fasta.out\" before removing CR1 annotated fragments.");log_file.write("Processing \"none_matched_seq.fasta.out\" before removing CR1 annotated fragments.")
# grep all lines containing LTR annotation and get list of containing fragment, then do the same for CR1
subprocess.call("grep \"LTR\" none_matched_seq.fasta.out | awk \'{print $5}\' | sort -n | uniq -c | awk \'{print $2}\' > LTR_matches.tmp", shell=True)
subprocess.call("grep \"CR1\" none_matched_seq.fasta.out | awk \'{print $5}\' | sort -n | uniq -c | awk \'{print $2}\' > CR1_matches.tmp", shell=True)
# cat these and extract only the unique lines, append this to LTR_matches (to account for unique CR1 lines) and extract the duplicated (LTR matches) only
subprocess.call("cat LTR_matches.tmp CR1_matches.tmp | sort -n | uniq -u >> LTR_matches.tmp", shell=True)
subprocess.call("sort -n LTR_matches.tmp | uniq -d > LTR_matches_only.tmp", shell=True)
ltr_annotation_only = open("LTR_matches_only.tmp").read().rstrip("\n").split("\n")

# process the LTR annotated positions, write to file, add to updated positions and sort positions
print("Fragments annotated as LTR retrotransposons identified. Writing to file.");log_file.write(" Fragments annotated as LTR retrotransposons identified. Writing to file. ")
ltrmatch = open("LTR_matches_to_add.tmp", "w")
if ((len(ltr_annotation_only) > 0) and ltr_annotation_only[0]):
	for hit in ltr_annotation_only:
        	if ((hit.split(":"))[0].count("_")) > 1:
                	chromo = str(hit.split("_")[0]) + "_" + str(hit.split("_")[1]) + "\t"
                	st = str((hit.split("_")[2]).split(":")[0]) + "\t"
                	end = str((hit.split("_")[2]).split(":")[1]) + "\t"
                	strand = str(hit.split("_")[3]) + "\tBR\n"
                	ltrmatch.write(chromo + st + end + strand)
        	else:
                	chromo = str(hit.split("_")[0]) + "\t"
                	st = str((hit.split("_")[1]).split(":")[0]) + "\t"
                	end = str((hit.split("_")[1]).split(":")[1]) + "\t"
                	strand = str(hit.split("_")[2]) + "\tBR\n"
                	ltrmatch.write(chromo + st + end + strand)
ltrmatch.close()
if os.path.isfile("validated_positions_update_01.txt"):
	subprocess.call("cat validated_positions_update_01.txt LTR_matches_to_add.tmp | sort -k1,1 -k2,2n > validated_positions_update_02.txt", shell=True)
else:
	subprocess.call("cat best_hits.txt LTR_matches_to_add.tmp | sort -k1,1 -k2,2n > validated_positions_update_02.txt", shell=True)

#process the CR1 fragments
print("Analysing fragments with LTR retrotransposon and CR1 annotation.");log_file.write("Analysing fragments with LTR retrotransposon and CR1 annotation.")
#extract all fragment names for LTRs and return all these sequences from file
subprocess.call("for i in `grep \"LTR\" none_matched_seq.fasta.out | awk \'{print $5}\' | sort -n | uniq -c | awk \'{print $2}\'`;do grep $i none_matched_seq.fasta.out >> all_ltr_annotations.tmp;done", shell=True)

if os.path.isfile("CR1_hits.tmp"):
	subprocess.call("rm CR1_hits.tmp", shell=True)
#if any of these sequences also have a CR1 annotation, retain these to analyse with the cr1_fragment_remover
if int(subprocess.check_output("grep -m 1 \"CR1\" all_ltr_annotations.tmp | wc -l", shell=True)) > 0:
        subprocess.call("for i in `grep \"CR1\" all_ltr_annotations.tmp | awk \'{print $5}\' | sort -n | uniq -c | awk \'{print $2}\'`;do grep $i all_ltr_annotations.tmp | awk \'{print $1\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$10\"\\t\"$11}\'>> CR1_hits.tmp;echo \"\" >> CR1_hits.tmp;done", shell=True)
        subprocess.call("cp CR1_hits.tmp nm_seq_CR1_annotations.txt", shell=True)
        static_functions.cr1_fragment_remover("nm_seq_CR1_annotations.txt", "validated_positions_update_02.txt", "03", log_file, "BR")
        subprocess.call("awk \'{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\tBR\"}\' validated_positions_update_03.txt > updated_positions_after_validation.txt", shell=True)
else:
        print("No LTR/CR1 annotations found.\n");log_file.write("No LTR/CR1 annotations found.\n")
        subprocess.call("awk \'{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\tBR\"}\' validated_positions_update_02.txt > updated_positions_after_validation.txt", shell=True)

static_functions.seq_extract("updated_sequences_after_validation.fasta", (static_functions.list_initial_formatter("updated_positions_after_validation.txt")), ref_gen_dict, "Y")
subprocess.call("rm -r *.fasta.*", shell=True)

end = subprocess.check_output("date | awk \'{print $3 \" \" $2 \" \" $6 \"  \"  $4}\'", shell=True)
subprocess.call("cp updated_sequences_after_validation.fasta ../refBLAST_final_identified_sequences.fasta", shell=True)
subprocess.call("cp updated_positions_after_validation.txt ../refBLAST_final_identified_positions.txt", shell=True)
os.chdir("../")
print("\nProcess complete.\nProcess began: " + str(startdate) + "Process finished: " + str(end) + "\n");log_file.write("\nProcess complete.\nProcess began: " + str(startdate) + "Process finished: " + str(end) + "\n")
log_file.close()