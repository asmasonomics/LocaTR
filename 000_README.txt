LocaTR README

Author: Andrew Mason
Date:	16th June 2016

Thank you for downloading LocaTR. Please cite the original paper in all work:

Mason AS, Fulton JE, Hocking PM & Burt DW (in press), A New Look at the LTR retrotransposon content of the 
	chicken genome, BMC Genomics


INSTRUCTIONS

LocaTR is designed to run (mostly) on a LINUX server. 

All scripts and directories within your download should be stored in one place, usually called "LocaTR" or similar.
You should then run the python script "000_modify_paths.py" from within this LocaTR directory.
This quickly goes through all script files and sets the path definitions to the download location.
You should then be able to run any script (given the appropriate path) successfully within the LINUX server. 

As will become clear, some scripts are for use on Windows machines. LTR_STRUC is a Windows executable, so must
be run in Windows. RetroTector can be run on Windows or Mac, but the preprocessing scripts I have written are
for Windows - some of the hard coded paths can be changed if you wish to use a Mac.

The LocaTR pipeline runs in 4 stages:
	1) Genome Preprocessing (detailed in 100_Genome-Processing_README.txt)
	2) Structurally Intact Element identification (detailed in 200_SIE-Identification_README.txt)
	3) Homology based identification (detailed in 300_Homology-Identification_README.txt)
	4) Validation (detailed in 400_Validation_README.txt)
	
Each of these README files will guide you through the appropriate methodologies. Stage 1 must be completed first.
Stages 2 and 3 can be run together. Stage 4 must be done following completion of both stages 2 and 3.


REQUIREMENTS

Most analysis and file processing occurs in a LINUX environment using BASH. Python scripts are most commonly used
and so Python 2.7 should be on the system. Scripts may run with Python 3 and above, but they have not been written
with Python 3 compatability in mind.

Individual program executables are not included within this release, but should be downloaded separately and installed.
In the first instance this should include Command line BLAST and HMMER version 3. Some helpful instructions are 
given for the installation of LTR_STRUC and RetroTector - the more complex programs. 


HELP

All documentation contained herein should allow you to run the LocaTR pipeline to identify LTR retrotransposons within
assembled genomes. However, there may be some compatability issues over time. If issues are specific to various 
identification programs please consult their manuals in the first instance. If it is due to the scripting or fomatting
of the output files, please feel free to contact Andrew Mason for help - andrew.mason@roslin.ed.ac.uk
He is also able to give tips on validation - the most difficult step.