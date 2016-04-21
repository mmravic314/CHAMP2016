## Marco Mravic DeGrado lab april 2016

# input 1: sequence string to thread onto an ideal helix
# input 2: path to 40 residue idealized alpha helix .pdb
# input 3: path to directory to put output/intermediate files
# input 4: path to rosetta main/

# output directory, input 3 path
# output 1: rosetta resfile for input sequence, starting chain A residue 1
# output 2: threaded and repacked helix, this repacking should't be taken seriously --just clearing clashes

### DEPENDENCIES: Rosetta compiled binaries, ProDy & all its dependencies (biopython, numpy)

# example command line
# python ~/CHAMP/bin/thread_repack.py PTGVIIGSIIAGILLLLALVAILWKL ~/bin/40_allAlaIdealHelix.pdb ~/CHAMP/a2/thread/ ~/rosetta/


import sys, os, subprocess as sp
from prody import *



## organize filepathes, inputs/outputs. load input pdb into prody

inSeq 	= sys.argv[1]
inPDB 	= parsePDB( sys.argv[2] )	

oDir 	= sys.argv[3]
if not os.path.exists( oDir ):
	os.mkdir( oDir )

rosiPath= sys.argv[4]
rpckPath=os.path.join( rosiPath, 'source/bin/fixbb.linuxgccrelease' )  

###

### 

# write resfile for sequence
txt 	= 'START\n'
resNum 	= 1
for res in inSeq:
	txt += '%d A PIKAA %s\n' % ( resNum, res )
	resNum += 1

outResF = os.path.join( oDir, 'resfile' )
oFile 	= open( outResF, 'w')
oFile.write( txt )
oFile.close()
# done writing res file

# trim and rewrite helix to length of input string
h_len	= len( inSeq ) 
if h_len > 40: 
	print "ERROR: input sequence longer than helix template input!!"

print h_len

cutPdb 	= inPDB.select( 'resnum 1 to %d' % h_len )

oPdb	= os.path.join( oDir, 'threadedHelix.pdb' )
writePDB( oPdb, cutPdb )
# done rewriting helix

# call rosetta to repack
cmd 	= [ rpckPath, 
'-s', 		oPdb, 
'-resfile', outResF, 
'-out:prefix', oDir,
'-nstruct', '1', 
'-ex1', '-ex2',
'-overwrite', 
'-minimize_sidechains'
]
print
print cmd
print 

sp.call(cmd)

#os.remove( oPdb )
