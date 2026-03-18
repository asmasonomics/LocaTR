## python 000_modify_paths.py [-h] 

import argparse
parser = argparse.ArgumentParser(
	description="modify_paths iterates through all scripts and changes base paths to match current directory system",
	epilog="Author: Andrew Mason; Release: 15/06/2016; Contact: andrew.mason@roslin.ed.ac.uk")

import os

new_wd = str(os.getcwd())
old_wd = "full_path_to_LocaTR"

for script in os.listdir("."):
	if (script.endswith(".py") or script.endswith(".sh")) and not script.startswith("000_modify"):
		t = open(script).read().rstrip("\n")
		out = open(script, "w")
		out.write(t.replace(old_wd,new_wd))
		out.close()