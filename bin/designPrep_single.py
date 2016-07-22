## Marco Mravic DeGrado Lab CHAMp design April 2016
#
# Prep an membrane-oriented protein  for rosetta design
#

# input 1: path to pdb file
# input 2: path to rosetta main
# input 3: aligned pdb file of the oriented template helix fro threading/ making resfile

# output: each pdb input file will have a design directory
# 			each dir will have a 'clean' rosetta-format input, a resfile, and span file from input conformation; also has output dir

# Example command line 
# python ~/CHAMP/bin/designPrep_single.py 4Helix_scrDSD_BAsm.pdb ~/rosetta/


import sys, os, subprocess as sp
from prody import *
from PDButil import UnNatAA

rosiBase 	= sys.argv[2]
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')


cmdSpan = [ rosiSpanGen, 
		'-database', rosiDB, 
		'-in:file:s', sys.argv[1]
		]		
print cmdSpan
sp.call( cmdSpan )
