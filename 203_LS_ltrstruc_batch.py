## 203_LS_ltrstruc_batch.py
## read the documentation before running this!!!

import shutil
import sys
import os
import re
import time
import win32com
import win32com.client

start_time = time.time()

ls_run_choice = raw_input("Please choose LTR_STRUC run sensitivity. Either: \n    Enter a number between 1 and 10 (1 being sensitive/slow and 10 being fastest)  or\n    Enter \"d\" for default.")

#processes the standard LTR_STRUC folder set up to make appropriate backups and remove the unnecessary "flist.txt"
if os.path.isfile(r"C:\LTR_STRUC\flist.txt"):
	os.remove(r"C:\LTR_STRUC\flist.txt")
shutil.copy(r"C:\LTR_STRUC\five_p_end.txt", r"C:\LTR_STRUC\COPY-five_p_end")
shutil.copy(r"C:\LTR_STRUC\log.txt", r"C:\LTR_STRUC\COPY-log")
shutil.copy(r"C:\LTR_STRUC\pbs.txt", r"C:\LTR_STRUC\COPY-pbs")
shutil.copy(r"C:\LTR_STRUC\rt.txt", r"C:\LTR_STRUC\COPY-rt")
os.mkdir(r"C:\LTR_STRUC\results")

#function for restoring LTR_STRUC folder to the standard set up following each iteration.
def folder_reset():
	from glob import glob
	for out_file in glob(r"C:\LTR_STRUC\*.txt"):
		os.unlink(out_file)
	shutil.copy(r"C:\LTR_STRUC\COPY-five_p_end", r"C:\LTR_STRUC\five_p_end.txt")
	shutil.copy(r"C:\LTR_STRUC\COPY-log", r"C:\LTR_STRUC\log.txt")
	shutil.copy(r"C:\LTR_STRUC\COPY-pbs", r"C:\LTR_STRUC\pbs.txt")
	shutil.copy(r"C:\LTR_STRUC\COPY-rt", r"C:\LTR_STRUC\rt.txt")
	return

import win32ui
def WindowExists(windowname):
		try:
			win32ui.FindWindow(None, windowname)

		except win32ui.error:
			return False
		else:
			return True
		
#iterates through all sequences in the given directory
for seq_file in os.listdir(r"C:\LTR_STRUC\sequences"):
	seq_name = (os.path.splitext(seq_file))[0]   #removes the extension from the file (LTR_STRUC will not run unless the seq file ends with ".txt")
	shutil.copy((r"C:\LTR_STRUC\sequences\\" + seq_file), (r"C:\LTR_STRUC\input\\" + seq_name + ".txt"))  #moves to input and adds ".txt"
	run_seq = open("flist.txt", "w")
	run_seq.write(seq_name + ".txt")
	run_seq.close()

	#LTR_STRUC executable is called and the defined keyboard entries are performed
	shell = win32com.client.Dispatch("WScript.shell")
	shell.Run(r"C:\LTR_STRUC\LTR_STRUC_1_1.exe")
	
	process = False
	while process == False:
		if WindowExists(r"C:\LTR_STRUC\LTR_STRUC_1_1.exe"):
			time.sleep(2)
			process = True
	
	shell.AppActivate(r"C:\LTR_STRUC\LTR_STRUC_1_1.exe")
	shell.SendKeys("y", 0)
	shell.SendKeys("{Enter}", 0)
	time.sleep(1)

	if ls_run_choice.isdigit():
		shell.SendKeys("n", 0)
		shell.SendKeys("{Enter}", 0)
		time.sleep(1)
		shell.SendKeys(str(ls_run_choice), 0)
		shell.SendKeys("{Enter}", 0)
		time.sleep(1)
		shell.SendKeys("{Enter}", 0)

	else:
		shell.SendKeys("y", 0)
		shell.SendKeys("{Enter}", 0)
		time.sleep(1)
		shell.SendKeys("{Enter}", 0)

	
	process = True
	while process == True:
		if WindowExists(r"C:\LTR_STRUC\LTR_STRUC_1_1.exe"):
			process = True
		else:
			process = False
			time.sleep(5)

	from glob import glob   #glob module performs name matching as would be expected in UNIX - allowing use of wildcards like *
	test = 0

	#counts all *_rprt.txt output files - if there aren't any this section will skip for a partciular input and the folder will reset ready for the next one
	for rprt_file in glob(r"C:\LTR_STRUC\*_rprt.txt"):
		test += 1 
	if test > 0:        
		j = 1
		for rprt_file in glob(r"C:\LTR_STRUC\*_rprt.txt"):
			content = open(rprt_file, "rU").read().rstrip("\n").split("\n")    #use of rU is Python 2* - treats the \r\n at the end of output file lines as \n

			#gets contig name
			name = seq_name.replace("_", "-")
			
			#finds the defining line for element length and extracts the value
			length = ""
			for line in content:
				if line.startswith("OVERALL LENGTH OF TRANSPOSON:"):
					length = line.replace("OVERALL LENGTH OF TRANSPOSON:", "").replace("bp", "").replace(" ", "")

			file_name = name + "_hit" + str(j).zfill(3) + "_" + str(length)   #creates an informative output style

			#this section writes the sequence to file following the format in the *_rprt.txt file
			#the file has the line shown below, then an empty line, then the sequence starts and runs on continuous lines until the next empty line
			#the while loop keeps writing until there is an empty line

			i = 0
			for line in content:
				out_file = open((file_name + ".fasta"), "a")
				if line.startswith("LENGTH OF PUTATIVE 5' LTR:"):
					q = i + 2
					prime5 = line.replace("LENGTH OF PUTATIVE 5' LTR:", "").replace("bp", "").replace(" ", "")
					prime3 = (content[q]).replace("LENGTH OF PUTATIVE 3' LTR:", "").replace("bp", "").replace(" ", "")
					header = ">" + str(prime5) + "-" + str(prime3)
					out_file.write(header + "\n")
				if line.startswith("COMPLETE SEQUENCE OF PUTATIVE TRANSPOSON:"):
					k = i + 2
					while (content[k]) != "":
						out_file.write(content[k])
						k+=1
				else:
					i+=1
				out_file.close()
			j+=1            
	
		for fasta_file in glob(r"C:\LTR_STRUC\*.fasta"):
			shutil.move(fasta_file, r"C:\LTR_STRUC\results\\")      #moves all .fasta files to a results directory for later use
	folder_reset()
	
#Following all iterations the LTR_STRUC folder is completely reset to the default setup - apart from the presence of the results folder!
os.remove(r"C:\LTR_STRUC\COPY-five_p_end")
os.remove(r"C:\LTR_STRUC\COPY-log")
os.remove(r"C:\LTR_STRUC\COPY-pbs")
os.remove(r"C:\LTR_STRUC\COPY-rt")
shutil.rmtree(r"C:\LTR_STRUC\input")
os.mkdir(r"C:\LTR_STRUC\input")
open(r"C:\LTR_STRUC\flist.txt", "a").close()

print("This process took " + "{0:.3f}".format(time.time() - start_time) + " seconds.")
