Homology Identification README file

Author: Andrew Mason
Date:	16th June 2016

All steps in 100_Genome-Processing_README.txt should have been completed before analysis begins.
Three Homology approaches will be used. Move into the directory 01_Homology/ for all analyses

1) RepeatMasker (RM)

RepeatMasker is the standard tool for identifying repetitive elements. For this type of annotation the search effort 
can be reduced by ignoring low complexity regions (-nolow flag) and by specifying a taxon area, even if it is still 
quite vague (-species vertebrates for instance). A good standard search would be:

	(i) repeatmasker -species vertebrates -nolow species_name_processed.fa &> /dev/null &

Positions can be extracted using a simple shell script, which also extracts the sequence:

	(ii) ./301_extract_RM_positions.sh species_name_processed.fa.out species_name_processed.fa

If you have a custom library of known elements not in RepBase then this can also be used. Be careful with this, as RM 
does not like long header names (must be less than 50 characters including ">") and this must include a repeat class 
definition. For us, "#LTR" is enough. So, a custom library must be a fasta file of sequences, with the names formatted 
to be: ">seq_name#LTR".  
The analysis can be run as so:

	(iii) repeatmasker -lib custom.fa -nolow species_name_processed.fa &> /dev/null &
	
Following both of these analyses you only need to keep the files ending .out and .tbl. Even the .tbl is not necessary, 
but it shows nice summary information about general repeat content etc.

Be careful with custom libraries, as any non-LTR elements in the custom library will find all of those. In mammalian 
genomes, the addition of a LINE would hugely increase annotated repeat content, but they are mostly going to be false 
positives in LTR retrotransposon identification. To ensure no non-LTR elements slip through, it is often prudent to run 
a second RM (using the -species vertebrate flag in this case) to check for non-LTR elements and remove them.

2) RefBLASTsearch

This stage is the most comprehensive for identifying specific ranges of elements for your species. Provided in this pipeline
is a fasta file of reference LTR retrotransposon sequences gathered from GyDB, augmented by sequences I found increased the
scope of this search - 303_refBLAST_full_ref_seq_from_GyDB_and_extras.fas. This file can be added to by the user in normal
fasta format. I find adding the sequence on one line only is always safest for cross-compatability.

From your homology directory run this command, specifying full paths:

	(v) 302_refBLASTsearch.py 303_refBLAST_full_ref_seq_from_GyDB_and_extras.fas species_name_processed.fa &> refBLAST.log &
	
	By default 8 threads are used for the BLAST run. This can be specified with the --threads command within this line.

This process creates its own results files within its own directory, and then copies out the important results files into the
original directory at the end. During the refBLASTsearch, identified elements are checked for annotation to other repeat classes
and incorrect matches are removed.

3) ReDoSt

ReDoSt detects DIRS1 elements within genomes using a database of existing elements contained within the installation. 
The program is an easy install, but needs Python 2.7 with NumPY and BioPython modules, and a standalone BLAST program. 
The genome to be studied should be put in the folder "genomes" in its own directory (use the same name, a softlink 
from the species_name_processed.fa to species_name/species_name.fasta is good - must use fasta end term). The program 
is then run (see below) and produces a directory in the results directory with results files. The file "annotCoGenom.tmp"
will be created if there are any hits. Sometimes there will not be, there has been significant lineage dependent sorting 
of these elements in Eukaryotes. 

	(vi) ReDoSt folder_where_genome_is > log.txt 

Any hits will have positions, sequences, and containing 10kb genomic fragments, for the DIRS element pol gene, but not 
the full element. This pipeline does not identify the full DIRS element, but simply extracts the annotated repeat region.
This will give an understatemenet on DIRS content. But these positions can be analysed further if the user wishes.
Whilst the script is used mainly to extract positions, it also does some basic validation using E-value threshold cut offs.

	(vii) 304_extract_dirs_positions.py -p PREFIX annotCoGenom.tmp
