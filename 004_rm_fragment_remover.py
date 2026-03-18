### python 004_rm_fragment_remover.py [-h] [--ori_seq_file ORI_SEQ_FILE] [-p] <rm_out> <ori_pos>

## define user arguments
from __future__ import division
import argparse
parser = argparse.ArgumentParser(
	description="rm_fragment_remover takes the output from RepeatMasker and removes fragments and sections with unwanted annotation.",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("rm_out", help="Specify RepeatMasker output file ending in .out")
parser.add_argument("ori_pos", help="Original positions file used to get sequences for RM run")
parser.add_argument("--ori_seq_file", help="Specify original fasta file pre editing for RM, if relevant (edit for too long seq headers for example)")
parser.add_argument("-p", "--prefix", action="store_true", help="if used program will use the first bit of the rm_out name (before first underscore) as prefix")
usr_args = parser.parse_args()

import os
import subprocess
import sys
sys.path.append("full_path_to_LocaTR/00_pipeline")
import static_functions

p=""
if usr_args.prefix:
        p = ((usr_args.rm_out).split("/")[-1]).split("_")[0] + "_"

# process RM out file if you had to edit seq header names before running RM analysis
if usr_args.ori_seq_file:
        i=1
        for head in (subprocess.check_output("grep \">\" " + usr_args.ori_seq_file, shell=True)).rstrip("\n").split("\n"):
                # use of \<\> forces exact match with sed, but need the single quotation marks
                subprocess.call(("sed -i \'s/\\<seq" + str(i) + "\\>/" + head.replace(">", "") + "/g\' " + usr_args.rm_out), shell=True)
                i+=1

# extract all hits which have annotations for CR1, LINE, DNA (transposon) or SINE, including all other annotations for that hit.
if os.path.isfile(p + "unwanted_fragment_hits.tmp"):
        os.remove(p + "unwanted_fragment_hits.tmp")
subprocess.call("for i in `egrep \"CR1|LINE|DNA|SINE\" " + usr_args.rm_out + " | awk \'{print $5}\' | sort -n | uniq -c | awk \'{print $2}\'`;do grep $i " + usr_args.rm_out + " | awk \'{print $1\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$10\"\\t\"$11}\'>> " + p + "unwanted_fragment_hits.tmp;echo \"\" >> " + p + "unwanted_fragment_hits.tmp;done", shell=True)

# process these hits with extra fragments and create a version of the original position file which lacks these hits using uniq
subprocess.call(("sed /^$/d " + p + "unwanted_fragment_hits.tmp | awk \'{print $2}\' | sort -n | uniq | sed s/:/\\\t/g | sed s/_/\\\t/g > " + p + "hits_with_non_ltr_fragments.txt"), shell=True)
subprocess.call(("cat " + p + "hits_with_non_ltr_fragments.txt " + usr_args.ori_pos + " | sort -n | uniq -u > " + p + "original_pos_no_non-ltr_annotations.txt"), shell=True)

# process the non-LTR annotations and extract sections of fragment that are LTR to be remerged with modified original pos list
raw_fragment_data = (open(p + "unwanted_fragment_hits.tmp").read().rstrip("\n")).split("\n\n")
fragment_data = []
# split data so each element of each individual fragment is available
for a in raw_fragment_data:
        x = []
        for b in (a.split("\n")):
                y = []
                for c in (b.split("\t")):
                        y += [c]
                if not ((y[5].startswith("Simple_repeat")) or (y[5].startswith("Satellite")) or (y[5].startswith("tRNA")) or (y[5].startswith("ncRNA")) or (y[5].startswith("rRNA")) or (y[5].startswith("Low")) or (y[5].startswith("Unknown"))): 
                        x += [y]
        fragment_data += [x]


new_hits = []
for hit in fragment_data:
        new_hit = ""
        if ((len(hit) > 1) and (str(hit).count("LTR") > 0)):
                x = len(hit) - 1
                #define key parts of fragment for rewriting
                if ((hit[0][1].split(":"))[0]).count("_") > 1:
                        chromo = str(((hit[0][1]).split("_"))[0]) + "_" + str(((hit[0][1]).split("_"))[1])
                        start_pos = (((hit[0][1]).split("_"))[2].split(":"))[0]
                        end_pos = (((hit[0][1]).split("_"))[2].split(":"))[1]
                        strand = ((hit[0][1]).split("_"))[3]
                else:
                        chromo = ((hit[0][1]).split("_"))[0]
                        start_pos = (((hit[0][1]).split("_"))[1].split(":"))[0]
                        end_pos = (((hit[0][1]).split("_"))[1].split(":"))[1]
                        strand = ((hit[0][1]).split("_"))[2]

                # get the indices for the LTR and non-LTR annotations to see if the search space can be narrowed
                ann = []
                for frag in hit:
                        ann += [str(frag[5])]
                ltr_match_indices = []
                non_ltr_match_indices = []
                for i,j in enumerate(ann):
                        if j.startswith("LTR"):
                                ltr_match_indices += [i]
                        else:
                                non_ltr_match_indices += [i]

                # see if two or more LTR annotations are adjacent or not - if not they might need to be split apart into two separate hits
                divided_ann = False
                if (len(ltr_match_indices) > 1):
                        k=0
                        while k < (len(ltr_match_indices) - 1):
                                if ((ltr_match_indices[k] + 1) < ltr_match_indices[k+1]):
                                        divided_ann = True
                                k+=1

                from itertools import groupby
                from operator import itemgetter
                range_list = []
                if divided_ann == False:        
                        # narrow processing space if there are runs of non LTR annotations at the start or end of a particular fragment
                        i_range = []
                        if (((ltr_match_indices[0] == 0) or (ltr_match_indices[0] == 1)) and ((ltr_match_indices[-1] == x) or (ltr_match_indices[-1] == x-1))):
                                i_range = range(0,x)
                        elif ((ltr_match_indices[0] >= 2) and ((ltr_match_indices[-1] == x) or (ltr_match_indices[-1] == x-1))):
                                i_range = range(((ltr_match_indices[0])-1),x)
                        elif (((ltr_match_indices[0] == 0) or (ltr_match_indices[0] == 1)) and (ltr_match_indices[-1] <= x-2)):
                                i_range = range(0,((ltr_match_indices[-1]) + 1))
                        else:
                                i_range = range(((ltr_match_indices[0])-1),((ltr_match_indices[-1]) + 1))
                        range_list += [i_range]
                else:
                        for k, g in groupby(enumerate(ltr_match_indices), lambda (i,x):i-x):
                                part_range = map(itemgetter(1), g)
                                # part range is not at the start of the indices and finishes before the end
                                if ((part_range[0] > 0) and (part_range[-1] < (len(ann) - 1))):
                                        range_list += [range((part_range[0] - 1),(part_range[-1] + 1))]
                                # part range starts at the start, and therefore must finish before the end (to be in this condition)
                                elif part_range[0] == 0:
                                        range_list += [range(0,(part_range[-1] + 1))]
                                # part range starts in the middle of the range and finishes at the end
                                elif (part_range[-1] == (len(ann) - 1)):
                                        range_list += [range((part_range[0] - 1),len(ann))]

                for i_range in range_list:
                        ch_hit = False
                        for i in i_range:
                                if ((hit[i][5]).startswith("LINE") or (hit[i][5]).startswith("DNA") or (hit[i][5]).startswith("SINE")):
                                        st = hit[i][2]
                                        end = hit[i][3]
                                        sw = hit[i][0]
        
                                        if ((hit[i+1][5]).startswith("LTR")):
                                                # non LTR element is before LTR element and overlaps, so compare SW scores
                                                if end > hit[i+1][2]:
                                                        if int(sw) > int(hit[i+1][0]):
                                                                new_start_pos = int(start_pos) + int(end)
                                                                new_hit = (str(chromo) + "\t" + str(new_start_pos) + "\t" + str(end_pos) + "\t" + strand + "\tBR")
                                                        else:
                                                                new_start_pos = int(start_pos) + int(hit[1][2])
                                                                new_hit = (str(chromo) + "\t" + str(new_start_pos) + "\t" + str(end_pos) + "\t" + strand + "\tBR")
                                                        ch_hit = True

                                                # non LTR element is before LTR element and doesn't overlap
                                                elif end < hit[i+1][2]:
                                                        new_start_pos = int(start_pos) + int(end)
                                                        new_hit = (str(chromo) + "\t" + str(new_start_pos) + "\t" + str(end_pos) + "\t" + strand + "\tBR")
                                                        ch_hit=True

                                        # non LTR element is after LTR element
                                        if ((i>0) and ((hit[i-1][5]).startswith("LTR"))):
                                                # element hasn't been modified already
                                                if ch_hit == False:
                                                        if st < hit[(i_range[-1])][3]: 
                                                                if int(sw) > int(hit[(i_range[-1])][0]):
                                                                        new_end_pos = int(start_pos) + (int(st) - 1)
                                                                        new_hit = (str(chromo) + "\t" + str(start_pos) + "\t" + str(new_end_pos) + "\t" + strand + "\tBR")
                                                                else:
                                                                        new_end_pos = int(start_pos) + int(hit[x-1][3])
                                                                        new_hit = (str(chromo) + "\t" + str(start_pos) + "\t" + str(new_end_pos) + "\t" + strand + "\tBR")
                                                        else:
                                                                new_end_pos = int(start_pos) + (int(st) - 1)
                                                                new_hit = (str(chromo) + "\t" + str(start_pos) + "\t" + str(new_end_pos) + "\t" + strand + "\tBR")
                                                else:
                                                        new = new_hit
                                                        if st < hit[(i_range[-1])][3]: 
                                                                if int(sw) > int(hit[(i_range[-1])][0]):
                                                                        new_end_pos = int(start_pos) + (int(st) - 1)
                                                                        new_hit = str(new.split("\t")[0]) + "\t" + str(new.split("\t")[1]) + "\t" + str(new_end_pos) + "\t" + strand + "\tBR" 
                                                                else:
                                                                        new_end_pos = int(start_pos) + int(hit[x-1][3])
                                                                        new_hit = str(new.split("\t")[0]) + "\t" + str(new.split("\t")[1]) + "\t" + str(new_end_pos) + "\t" + strand + "\tBR"
                                                        else:
                                                                new_end_pos = int(start_pos) + (int(st) - 1)
                                                                new_hit = str(new.split("\t")[0]) + "\t" + str(new.split("\t")[1]) + "\t" + str(new_end_pos) + "\t" + strand + "\tBR"

                new_hits += [new_hit]

        
#write newly edited hits to tmp file
added_hits = open(p + "fragments_removed.tmp", "w")
for hit in filter(None, new_hits):
        added_hits.write(hit + "\n")
added_hits.close()

subprocess.call("cat " + p + "fragments_removed.tmp " + p + "original_pos_no_non-ltr_annotations.txt | awk '{if ($3<$2) {print $1 \"\t\" $3 \"\t\" $2 \"\t\" $4 \"\t\" $5} else {print $0}}' > " + p + "fragments_and_ori.tmp", shell=True)
subprocess.call("sort -k1,1 -k2,2n " + p + "fragments_and_ori.tmp | uniq > " + p + "fragments_removed.txt", shell=True)
subprocess.call("rm " + p + "fragments_removed* " + p + "hits*", shell=True)