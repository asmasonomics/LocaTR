SIE Identification README file

Author: Andrew Mason
Date:	16th June 2016

All steps in 100_Genome-Processing_README.txt should have been completed before analysis begins.

Four structural identification programs will be used. Move into the 02_SIE/ directory and, for safety, create separate directories
for each identification methodology.

Whilst we recommend the use of LTR_STRUC, LTR Harvest, MGEScan_LTR and RetroTector, other LTR retrotransposon or repeat identification
software can be easily incorporated into the LocaTR pipeline. This may be necessary as we recommend a multi-tool approach and some 
programs are limited to certain architecture (LTR_STRUC to Windows, and RetroTector to Mac or Windows, not Linux). 
Results from any additional programs can be added as part of the validation stage using tab spaced positions files in the format:
Contig	Start_position	End_position	Strand	Program_code
Program code is any string. For example, the program code for LTR Harvest is LH. 
The format is similar to a modified BED, but positions are 1-indexed. Ensure contig names match those in the preprocessed genome.

Requirements:
LTR Harvest is within the Genome Tools package downloaded in the genome preprocessing steps. 
LTR_STRUC can be requested from the McDonald Lab at Georgia Tech - see http://www.mcdonaldlab.biology.gatech.edu/ltr_struc.htm for 
up to date contact information.
MGEScan_LTR is available from - https://sourceforge.net/projects/mgescan/
RetroTector can be requested from the Blomberg Lab at Uppsala University - see 
http://www.mybiosoftware.com/retrotector-retrotectorviewer-1-0-1-identification-characterization-proviral-sequences-vertebrate-genomes.html
for up to date contact information

If there are any problems obtaining these last three, please contact Andrew Mason - andrew.mason@roslin.ed.ac.uk

Program specific requirements are detailed within each section. 

1) LTR Harvest (LH)

Part of the Genome Tools package, this is powerful and very quick following construction of the suffix array in Genome Processing. 
The three options below for -minlenltr, -maxlenltr and -similar are based on my work. These do not have to be followed exactly, and 
the values may have very different effects dependent on the genome. Please refer to the latest LTR Harvest manual for other options. 
Program should be run, then the gff3 sorted.

	(i) gt ltrharvest -index species_name_processed.fa -minlenltr 80 -maxlenltr 2000 -similar 75 -gff3 species_name_LH_out.gff3
	(ii) gt gff3 -sort species_name_LH_out.gff3 > species_name_LH_out.sorted.gff3

The gff3 file keeps information on features (mostly just the LTR pair), but for validation we need the sequences and positions using 
the external positions. These data can be extracted using "201_extract_LH_positions.py"

	(iii) python 201_extract_LH_positions.py species_name_LH_out.sorted.gff3 species_name_processed.fa

Positions are unstranded at this point, so strandedness will be added later. All extracted sequences come from the +ve strand.
LTR Harvest creates a lot of false positives. Positions should be thoroughly validated before being used (see 400_Validation_README.txt)


2) LTR_STRUC (LS)

It is quite a basic piece of software and has issues with memory allocation and will not reverse complement sequences to search in both 
the forward and reverse orientations. It will also not perform batch processing. Here we use three scripts to largely automate the process.
The main LTR_STRUC program is a windows executable.

"Extras" you will need:
- win32com <Python-version> Python module 
- Windows Python - available at https://www.python.org/download/

a) Sequence pre-processing
LTR_STRUC has memory allocation problems when files are too large. A simple solution is to split your multi-sequence FASTA file up into 
chromosomes or smaller contigs. Creation of reverse complement sequences is also required  for LTR_STRUC. Create a directory called 
LS_seq and run:

	(iv) python 202_LS_seq_formatter.py [-h] input_seq.fasta

All generated sequences should be moved into a directory called "sequences" and copied to your Windows desktop ready for use in stage 2.

b) Executing LTR_STRUC
Copy the "sequences" directory from above into the LTR_STRUC directory. Copy the "203_LS_ltrstruc_batch.py" file into the LTR_STRUC 
directory. Then, move this whole directory onto the C:\ drive so that it executes from there. The LTR_STRUC directory must be at the 
position: C:\LTR_STRUC

Navigate into this directory and run 203_LS_ltrstruc_batch.py. The script will immediately ask for user input (just once, at the start) 
asking you to choose run sensitivity. This can either be a number in the range 1-10, where 1 is the most sensitive (also slowest) and 10 
is the quickest, or by entering the string "d" or "D" to represent "default". This is roughly equivalent to selecting 5 on the sensitivity 
scale. Press enter after you have typed your choice. If you choose an invalid character the program will display an error message and exit.

Your machine must remain on at all times during this stage. There is no ability to put this process in the background.
This will run for some time (large genomes may take 3-4 days depending on your machine) and will produce a directory called "results". Extract 
this and return it to the LINUX server (preferably) for the final step.

WARNINGS!
i) This process runs in the foreground. If you lock your PC the current sequence analysis will complete, but the next will not initiate. 
ii) As each new execution of LTR_STRUC_1_1.exe begins, the new window is in the foreground. As a result you can accidentally type in the window 
causing the analysis for that particular sequence to crash. This is annoying I know, but you don't have  to restart the whole analysis. Just 
check your output and see which sequence has been skipped and rerun the python script with just that sequence as input.

c) Creating a positions list
Change into results directory and run:

	(v) python 204_LS_ltrstruc_batch_pos_extract.py [-h] ref_genome

Following position extraction the positions are ordered and merged (in case there  are duplicates following detection of both orientations of 
the sequence(s)). This final positions file is a tab delimited file with columns: chromosome-start-end-strand-source
In this file the source is "LS" for LTR_STRUC. This is helpful if you later use further identification methods and merge all datasets. 


3) MGEScan_LTR (MGS)

Create a results directory "MGS_results", and within this directory make another called "genome". Enter the 
genome directory and process the reference genome into single fasta sequences using the script below:

	(vi) python 205_MGS_seq_formatter.py species_name_processed.fa
	
Following this the MGEScan_LTR program can be run:

	(vii) /path/to/run_MGEScan.pl -genome=/full/path/MGS_results/genome/ -data=/full/path/MGS_results/ -hmmerv=3 -program=L &> species_name_MGS.log &
	
Positions can be extracted from the ltr/ltr.out file using a simple positions script. This also extracts the sequences.

	(viii) ./206_extract_MGS_positions.sh MGS_ltr.out species_name_processed.fa


4) RetroTector (ReTe/RT)

Create a directory called RT_results and run the following script to prepare the ReTe sequences:

	(ix) python 207_rete_input_fasta_formatter.py species_name_processed.fa
	
Then move the sequences to the desktop of a Windows or Mac (the instructions below will be for Windows machine). The next bit is quite complex, so take your time!





