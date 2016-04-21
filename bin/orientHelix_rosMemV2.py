# Marco Mravic DeGrado Lab 2016

## Input 1: pdb of ideal helix threaded with target sequence, e.g. a5-integrin TM region, and oriented originally by OPM/PPM or EZ
## Input 2: path 2 Octopus txt file
## Input 3: path to rosetta main
## input 4: path to Rosetta scripts XML with protocol 

# puts all output and intermediate files to the directory in the path of the input helix

## Example command line
# 
# python ~/CHAMP/bin/orientHelix_rosMemV2.py ~/CHAMP/a5/RosMem_Target-input/alignedInput.pdb  ~/rosetta/ ~/CHAMP/bin/helix_Relax.xml ~/CHAMP/a5/thread/resfile 

# example sequence for OCTOPUS (http://octopus.cbr.su.se/), a5-integrin
# PLWIIILAILFGLLLLGLLIYILYKLG


import sys, os, subprocess as sp
from prody import *

#### I/Oa5`
inPDB		= sys.argv[1]


rosiBase 	= sys.argv[2]
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiOctPrl 	= os.path.join( rosiBase, 'source/src/apps/public/membrane_abinitio/octopus2span.pl')
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')

protocolPth = sys.argv[3]
resfile_path= sys.argv[4]

oDir 		= os.path.dirname( inPDB ) + '/'
####

## Clean DUm atoms away 
#inPDB_clean = inPDB[:-4] + '_clean.pdb'
#pdb2 		= parsePDB( inPDB ).select( 'protein' ).copy()
#writePDB( inPDB_clean, pdb2 )
##


## Make span file from input PDB

cmdSpan = [ rosiSpanGen, 
'-database', rosiDB, 
'-in:file:s', inPDB
]

sp.call( cmdSpan )
#sys.exit()

spanF = inPDB[:-4] + '.span'

### Relax

cmd = [  rosiScrps, 
'-parser:protocol', protocolPth, 			# Path to Rosetta script (see above)
'-in:file:s', inPDB,							# Input PDB structure
'-nstruct', '1', 							# Generate 1000 models
'-mp:setup:spanfiles', spanF,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix', oDir,
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

print 
print cmd
print 

sp.call( cmd )

sys.exit()

###

############# Viewing option below

## View in membrane

rosiMemPyMol = os.path.join( rosiBase, 'source/bin/view_membrane_protein.linuxgccrelease')
rosiMemPyMoX = os.path.join( rosiBase, 'source/bin/view_membrane_protein.linuxgccrelease')
pythonPymol	 = os.path.join( rosiBase, 'source/src/python/bindings/PyMOLPyRosettaServer.py')
#run /home/xray/rosetta/source/src/python/bindings/PyMOLPyRosettaServer.py

# python ~/CHAMP/bin/orientHelix_rosMem.py ~/CHAMP/a5/RosMem_Target-input/alignedInput.pdb ~/CHAMP/a5/RosMem_Target-input/a5-Octopus.txt ~/rosetta/ ~/CHAMP/bin/helix_Relax.xml 



'''
/home/xray/rosetta/source/bin/view_membrane_protein.linuxgccrelease -database /home/xray/rosetta/database/ -in:file:s example_inputs/1c3w_tr.pdb -mp:setup:spanfiles example_inputs/1c3w_tr.span -show_simulation_in_pymol 0 -keep_pymol_simulation_history 1



'''

print pythonPymol

cmdOG = [ rosiMemPyMol, 	'-database', rosiDB,
	'-in:file:s', inPDB,
	'-mp:setup:spanfiles', oSpan,
	'-show_simulation_in_pymol', '0',
	'-keep_pymol_simulation_history', '1', 
	'-out:overwrite'
	]

sp.call( cmdOG )



sys.exit()
#
#sp.call( rosiMemPyMol ) 

# extract membrane position from pdb output
for i in os.listdir( oDir ):
	if '_' not in i: continue

	path = os.path.join( oDir, i )
	## Clean DUm atoms away 
	inPDB_clean = inPDB[:-4] + '_clean.pdb'
	pdb2 		= parsePDB( inPDB ).select( 'protein' ).copy()
	writePDB( inPDB_clean, pdb2 )


	cmd = [ rosiMemPyMol, 	'-database', rosiDB,
	'-in:file:s', inPDB_clean,
	'-mp:setup:spanfiles', oSpan,
	'-show_simulation_in_pymol', '0',
	'-keep_pymol_simulation_history', '1'

	]

	print
	print ' '.join( cmd )
	print







##

'''  OPTIONS FROM 2015 PLOS paper
-parser:protocol membrane_relax.xml # Path to Rosetta script (see above)
-in:file:s 3PXO_tr_native.pdb # Input PDB structure
-nstruct 1000 # Generate 1000 models
-mp:setup:spanfiles 3PX0.span # Input spanfile
-mp:scoring:hbond # Turn on membrane hydrogen bonding
-relax:jump_move true # Allow jumps to move during relax
-packing:pack_missing_sidechains 0 # Wait to pack sidechains 
'''