#! /usr/bin/python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='input_blast', dest='i',
					type=str, required=True,
					help='file containing Blast results to be used as input')
#parser.add_argument('-c', '--coverage', metavar='coverage_set', dest='c',
					#type=float, required=True,
					#help='set coverage')
parser.add_argument('-l', '--plasmid_length', metavar='plasmid length', dest='l',
					type=int, required=True,
					help='the length of plasmid for blast')

args = parser.parse_args()

#set List according to plasmid length
plasmid_length = []
n=1
while n < args.l + 1:
	plasmid_length.append(0)
	n += 1
#print(len(plasmid_length))

INPUT_FILE = args.i
for line in open(INPUT_FILE, "r"):
	line.rstrip()
	info = line.split()
	if float(info[2]) > 0:
		if int(info[8]) < int(info[9]):
			for i in range(int(info[8]), int(info[9])+1):
				plasmid_length[i-1] = 1
		else:
			#print(info[8])
			for i in range(int(info[9]), int(info[8])+1):
				plasmid_length[i-1] = 1

#print(plasmid_length)
total = 0
ele = 0
while(ele < len(plasmid_length)):
	total = total + int(plasmid_length[ele])
	ele += 1

coverage = int(total)/args.l
print(INPUT_FILE, "\t", coverage)
#print('\n') 
