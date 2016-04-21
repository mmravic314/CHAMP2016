# Marco Mravic DeGrado Lab April 2016
# CHAMP design

# input 1: integrin input structure output from PPM/OPM server, pdb 2knr 
# input 2: structure of threaded pdb, aka the template integrin, e.g. a5
# input 3: output path directory for aligned pdb


# python ~/CHAMP/bin/alignTMs.py ~/CHAMP/Misc/2knc.pdb ~/CHAMP/a5/thread/threadedHelix_0001IN.pdb ~/CHAMP/a5/RosMem_Target-input/



## 'Known' TM region given the OPM/PPm prediction for the aii-b1 NMR structure
aii_TM_seq	= 'IWWVLVGVLGGLLLLTILVLAMWKVGFFK'
##


import sys, os, re
from PDButil import UnNatAA
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from prody import *


## Alignment info
matrix = matlist.blosum62
gap_open = -20
gap_extend = -0.1
##

oPath = os.path.join( sys.argv[3], 'alignedInput.pdb' )

# load and renumber aii-b chain
aii_strct	= parsePDB( sys.argv[1], chain='A', subset='bb' )
step = 1
for resi in aii_strct.iterResidues():
	resi.setResnum( step )
	step += 1 

aii_seq		= ''.join( [ UnNatAA[r] for r in aii_strct.select( 'ca' ).getResnames() ] )

in_strct 	= parsePDB( sys.argv[2] )
in_seq		= ''.join( [ UnNatAA[r] for r in in_strct.select( 'ca' ).getResnames() ] )


## Do alignment of aii structure (seq) and a5 structure (seq)
alns = pairwise2.align.globalds( aii_TM_seq, in_seq, matrix, gap_open, gap_extend)
top_aln = alns[0]
aln_aii, aln_in, score, begin, end = top_aln

print '\nalign ment found: a-IIb TM region to input helix'
print aln_aii+'\n'+aln_in
print 
## Grab index of starting TM region for input
start 		= re.search( r'-+\w', aln_aii ).end() 
end 		= re.search( r'\w+-', aln_aii ).end() -1

## selection string and Atom Group for aligned TM region of in-helix
in_selStr 	= 'ca resnum %d to %d' % ( start, end )
in_mobile 	= in_strct.select( in_selStr ).copy()

## selection string for aligned TM region of aiib helix
tmpMtch		= re.search( r'%s' % aii_TM_seq, aii_seq )
aii_srt, aii_end = tmpMtch.start() + 1, tmpMtch.end()
aii_selStr	= 'ca resnum %d to %d' % ( aii_srt, aii_end ) 
aii_target 	= aii_strct.select( aii_selStr ).copy()

# Align mobile in_helix to aiib- target helix based on this sequence alignment
transMat	= superpose( in_mobile, aii_target )[1]
#print 'CA rmsd over %d atoms' % len( aii_TM_seq ), calcRMSD( in_mobile, aii_target )
outPDB 		= applyTransformation( transMat, in_strct ) 
print 'wrote aligned helix'

writePDB( oPath, outPDB )




