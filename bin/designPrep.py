## Marco Mravic DeGrado Lab CHAMp design April 2016
#
# Prep a directory of backbone models to design. 
#

# input 1: path to directory containing 'match_x.pdb files'
# input 2: path to rosetta main
# input 3: aligned pdb file of the oriented template helix fro threading/ making resfile

# output: each pdb input file will have a design directory
# 			each dir will have a 'clean' rosetta-format input, a resfile, and span file from input conformation; also has output dir

# Example command line 
# python ~/CHAMP/bin/designPrep.py ~/CHAMP/a5/bb2Design/ ~/rosetta/

import sys, os, subprocess as sp
from prody import *
from PDButil import UnNatAA

rosiBase 	= sys.argv[2]
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')


for f in os.listdir( sys.argv[1] ):
	print f
	if f[-4:] != '.pdb' or 'membrane' in f: continue
	
	## i/o for each input match
	inP 		= os.path.join( sys.argv[1], f)
	inPDB 		= parsePDB( inP )
	print 'entering', inP

	oDir		= inP[:-4]
	oPDBpath	= os.path.join( oDir, f )
	oPResF		= os.path.join( oDir, 'resfile' )
	
	if not os.path.exists( oDir ):
		os.mkdir( oDir ) 

	##

	# Write res file, using integrin input sequence on chain A, and anything on chain X
	txt = 'start\n* X ALLAAxc\n'
	for r in inPDB.iterResidues():

		ch, num, resi = r.getChid(), r.getResnum(), r.getResname()
		if ch == 'A':
			txt += '%d %s PIKAA %s\n' % ( num, ch, UnNatAA[resi] )
	resF = open( oPResF, 'w')
	resF.write( txt )
	resF.close()

	sp.call( ['mv', inP, oDir ] )

	cmdSpan = [ rosiSpanGen, 
		'-database', rosiDB, 
		'-in:file:s', oPDBpath
		]		
	print cmdSpan
	sp.call( cmdSpan )
