## python rete_positions_formatter.py [-h] [--prefix PREFIX] <SQL_username> <SQL_password> <SQL_host> <SQL_port> <SQL_database>
## must be run from server with SQL access and with the environmental variable $DBUSER defined as login details for your particular server

import argparse

parser = argparse.ArgumentParser(
	description="rete_positions_formatter accesses the ReTe output from the SQL database and processes it to extract putative element positions and sequences.",
	epilog="Author: Andrew Mason; Release: 03/11/14; Contact: andrew.mason@roslin.ed.ac.uk")

parser.add_argument("SQL_username", help="Specify SQL username for accessing ReTe output")
parser.add_argument("SQL_password", help="Specify SQL password for accessing ReTe output")
parser.add_argument("SQL_host", default="localhost", help="Specify SQL host for accessing ReTe output (default = \"localhost\")")
parser.add_argument("SQL_port", type=int, default=0, help="Specify SQL port for accessing ReTe output (default = 0)")
parser.add_argument("SQL_database", help="Specify name of SQL database which stores the ReTe output")
parser.add_argument("--prefix", default="RT", help="Specify a desired output file prefix (default = \"RT\")")
usr_args = parser.parse_args()

mysql_info = "mysql -h " + usr_args.SQL_host + " -P " + str(usr_args.SQL_port) + " -u " + usr_args.SQL_username + " -p" + usr_args.SQL_password + " " + usr_args.SQL_database + " -BN -e"
pre = usr_args.prefix

import sys
import re
import subprocess


### PART 1 - format positions file

ltr_info = (subprocess.check_output(mysql_info + " \"select id, FirstPos, LastPos from ltrs\" | awk \'{if($3>$2){print $1\" \t\"(($3+1)-$2)} else {print$1\" \t\"(($2+1)-$3)}}\'", shell=True)).rstrip("\n").replace("\r", "").split("\n")
ltr_data = {}
for line in ltr_info:
	ltr_data[(line.split("\t"))[0]] = (line.split("\t"))[1]

raw_data = subprocess.check_output(mysql_info + " \"select Chromosome, Start, End, Strand, LTR5id, LTR3id from chains\" | sed s/.fasta//g | awk \'{if($5~\"NULL\" && $6~\"NULL\") {print $0 \"\tRT\"} else {print $0\"\tRT_LTR\"}}\' | awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$7\"\t\"$5\"\t\"$6}\'", shell=True)
positions_split = (read_file.rstrip("\n").replace("\t", ", ")).split("\n")
x = []
for element in positions_split:
	x = (element.split(", "))
	positions = []
	if int(x[1]) > int(x[2]):
		st, end = 1, 2
		x[end], x[st] = x[st], x[end]
	if x[3] == "S":
		st, end = 5, 6
		x[end], x[st] = x[st], x[end]
	i=1
	for value in x:
		if i in [2, 3]:
			positions += [(int(value))]
		elif i == 4:
			if value == "P":
				positions += ["+"]
			else:
				positions += ["-"]
		elif i in [6, 7]:
			match = 0
			for key,val in ltr_data.iteritems():
				if key.startswith(value + " "):
					positions += [(int(val))]
					match = 1
			if match == 0:
				positions += [0]
		else:
			positions += [value]
		i+=1
	x += [positions]
x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))


### PART 2 - remove dubious padded chromosome hits

i=0
j=0
to_remove = []
for position in x:
	pad = r"_15k-padded_RL-"
	if re.search(pad, position[0]):
		max_len = (((position[0]).split("-"))[2]).replace(".", "").replace("kb", "").lstrip("0")
		
		if ((int(position[1]) - 15000) < 0) or ((int(position[2]) - 15000) > int(max_len)):
			to_remove += [i]
			j += 1
		else:
			x += [[((((position[0]).split("-"))[0]).replace("_15k", "")),(int(position[1]) - 15000),(int(position[2]) - 15000),(position[3]),(position[4]),(position[5]),(position[6])]]
			to_remove += [i]
	i+=1

print(str(j) + " padded elements were discarded due to padding overlap.")

for offset, index in enumerate(to_remove):
	index -= offset
	del x[index]



### PART 3 - merge positions and generate final list

i = 0
sort = None
x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))

while i < len(x) - 1:
	# performs a fresh sort at the start of iteration if the list order has changed
	if sort == True:
		x.sort(lambda l1, l2: cmp(l1[0], l2[0]) or (l1[0] == l2[0] and (cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))) or (l1[2] == l2[2] and cmp(l1[3], l2[3])))))
	
	if (x[i][0] == x[i + 1][0]) and (x[i][1] == x[i + 1][1]) and (x[i][2] == x[i + 1][2]):
		new = ((str(x[i])).replace("\'", "").replace("[", "").replace("]", "").replace("+", "m").replace("-", "m")).split(", ")
		new_pos = []
		count = 0
		for value in new:
			count += 1
			if count in [2, 3, 6, 7]:
				new_pos += [int(value)]
			else:
				new_pos += [value]
		x += [new_pos]
		x.pop(i+1)
		x.pop(i)
		sort=True
	else:
		i += 1
		sort = False 

outname = str(pre) + "_positions.txt"
outfile = open(outname, "w")
for pos in x:
	chromo = str(pos[0])
	if chromo.startswith("JH") or chromo.startswith("AADN"):
		chromo = chromo.replace("_", ".")
	outfile.write(chromo + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + str(pos[3]) + "\t" + str(pos[4]) + "\t" + str(pos[5]) + "\t" + str(pos[6]) + "\n")
outfile.close()

print("Succesfully created " + outname)
print("Total putative elements: " + str(len(x)))
tot_length = str(subprocess.check_output("awk \'{sum+=($3-$2)}END{print sum / 1000000, \"Mb\"}\'" + outname, shell=True))
print("Total sequence length: " + tot_length)