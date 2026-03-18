Validation README file

Author: Andrew Mason
Date:	16th June 2016

Elements generated from the 200_SIE-Identification_README.txt and 300_Homology-Identification_README.txt should be used following this README file.

SIEs are validated using Dfam nucleotide pHMMs, GyDB2.0 protein pHMMs and the tRNA pHMMs generated when running 102_extract_tRNA_seq.py
The user must also make some judgements here, as there is no hard coded cut off with many of the assignments. 
Given the information you must choose which intact elements to take forward. I provide some tips on making these calls below. 
Also, when you later merge positions files of multiple identification techniques you can choose to do further 
removal of putative elements by only taking those through that were identified in multiple analyses. 

1). VALIDATION OF SIEs

For each of the structural identification methods used, you need to individually validate the elements. Do this within the 02_SIE/
directory, but make separate validation directories - e.g. "LH_validation/" for each program.

You could combine lists at this point, but it can often generate huge files and you lose information on program-specific false positive rates etc.
You should make a new directory for each analysis, move into it and run the following command, using correct path definitions.

	(i) python 403_validate_SIEs.py -d 401_Dfam_nucl_pHMMs/ -p 402_GyDb_protein_pHMMs/ -t 02_SIE/tRNA_pHMMs/ -op OUT_PREFIX species_name_method_positions.txt species_name_method_seq.fa
	
Whilst running, this prints progress to stdout. On completion it generates three output files (all other files are removed - hence doing it in separate 
directories). The key file is OUT_PREFIX_potentially_validated_output.txt

This file can then be transferred into an excel spreadsheet for easier management. This could be done within linux but I find it is quicker, especially 
as files never usually get that huge. Select everything from row 2 downwards and sort, using column headers, by "Validated", then "Chromosome", then "Start"

From this point it is "just" a case of selecting the putative elements with sufficient validation.
The validation column shows if there is support from: 
	- Dfam with the "D" nomenclature, with the immediately following number showing where the first LTR annotation is within the list (sorted by E-value)
	- Protein analysis with the "P" nomenclature, where there is at least one protein domain match
	- RepeatMasker analysis with "R" nomenclature, where an element is given "R" if LTR annotation has the greatest element coverage
	- tRNA presence with "t" nomenclature, where an element is given "t" if a tRNA site is identified (lowest importance in validation)
Any elements marked "XX" can be immediately discarded.
	
The best scenario for any element is "D1+P+R", and these are usually validated at high confidence.
Watch out for LINE defintions within the RM. Is the proportion similar to the LTR proportion? If so, are there a low number of DNA motif matches?
Protein matches will ideally have a low E value score for RT, and the presence of other domains likely with higher E-values. 
The key here is to avoid bringing through LINEs due to their RT homology with LTR retrotransposons.

With P, or P+R, make a judgement on how low the protein scores are. These definitions often mean the pol gene has been disrupted/degraded.

With R alone, check the LTR coverage. Often better to ignore these. High coverage elements will be identified during the homology search. Low coverage ones,
or those with high LINE homology gives no element validation certainty.

t alone, just discard.

Lower D numbers (with or without protein or repeatmasker support), take care with these. When on their own (e.g. "D4") look at how many LTR motifs there are 
compared to other Dfam annotations. Likely that these are best to discard.

Highlight elements you wish to discard, and then copy columns B, C, D, E and G of those you want to keep into a notepad++ file. It is often a good idea to check
through the strand information at this point. Protein classifications will give the correct strand for the element for example.
Check to remove any special characters (\r mostly) and then save to file and transfer back to the server.

Once completed for all SIE programs used, transfer your validated positions files to the 03_Validated/ directory


2). SECONDARY BLAST

Once all validated positions lists are back on the server and within the 03_Validated/ directory, run the merger protocol:

	(ii) ./404_merge_positions_list.sh FILE1 FILE2 FILE3 .... FILEn

This will generate a single file of positions, where the fifth column shows the program(s) where the element was identifed. Change the name of the output file 
to something memorable.

You can manually select on elements with multiple program identifications at this stage if you wish. I do not suggest this, but it would give your (largely) 
reduced list of elements more confidence.

The next step is a secondary BLAST of the genome using just elements identified by the struc methods but not in the homology search. Make a new directory and
run 405_secondary_BLAST_analysis.py:

	(iii) python 405_secondary_BLAST_analysis.py -p OUT_PREFIX merged_SIE_positions.txt merged_homology_positions.txt /full/path/species_name_processed.fa
	
Note - you must use the full path to the processed reference genome, not a relative path.

As with the refBLAST in the homology searches, this will generate a list of positions. These can be merged with existing lists to generate a full list:

	(iv) ./404_merge_positions_list.sh [homology_program_pos_lists.txt] SIE_pos.txt secondary_BLAST_pos.txt


3). LAST STEPS

You now have all your positions! These can be extracted using:

	(v) python 001_seq_extract.py -u --prefix OUT_PREFIX merged_positions_file.txt species_name_processed.fa
	
You can then revert the positions and sequence contig/chromosome names back to their originals using this script for each file indivdiually (pos, then seq etc).
It can also be run on any of the intermediate files throughout the analysis if you want to check things.

	(vi) python 406_convert_back_to_original_contig_names.py -p OUT_PREFIX file_to_convert species_name_information.txt

Note - the "species_name_information.txt" was one of the files generated at the very start by the 101_format_genome_file.py script