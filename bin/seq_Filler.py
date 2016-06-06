# read a asta file of sequences with variable positions 'X' and return sequence with semi-random hydrophobics in that position
# input 1: fasta file
# output to standard out, pipe to new file

import sys
from numpy.random import random_integers


choices = ['L', 'L', 'L', 'L', 'L', 'L', 'A', 'F', 'V', 'I']

#v = random_integers(0, 9, 1)

with open( sys.argv[1] ) as fin:
	for i in fin:
		if i[0] == '>': 
			print i 
			continue

		print i
		seq = ''
		for r in i.rstrip():
			if r == '-': continue
			if r == 'X': 
				seq += choices[ random_integers(0, 9, 1) ]
				continue
			seq += r
		print seq
		print
		
