# Marco Mravic DeGrado Lab 2016

## Input 1: pdb of ideal helix threaded with target sequence, e.g. a5-integrin TM region, and oriented originally by OPM/PPM or EZ
## Input 2: path 2 Octopus txt file
## Input 3: path to rosetta main
## input 4: path to Rosetta scripts XML with protocol 

# puts all output and intermediate files to the directory in the path of the input helix

## Example command line
# 
# python ~/CHAMP/bin/orientHelix_rosMem.py ~/CHAMP/a5/RosMem_Target-input/orientedOPMA5.pdb ~/CHAMP/a5/RosMem_Target-input/a5-Octopus.txt ~/rosetta/ ~/CHAMP/bin/helix_Relax.xml 


# example sequence for OCTOPUS (http://octopus.cbr.su.se/), a5-integrin
# PLWIIILAILFGLLLLGLLIYILYKLG


import sys, os, subprocess as sp
from prody import *

#### I/Oa5`
inPDB		= sys.argv[1]
octF 		= sys.argv[2]

rosiBase 	= sys.argv[3]
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiOctPrl 	= os.path.join( rosiBase, 'source/src/apps/public/membrane_abinitio/octopus2span.pl')

protocolPth = sys.argv[4]

oDir 		= os.path.dirname( inPDB ) + '/'
####

## Clean DUm atoms away 
inPDB_clean = inPDB[:-4] + '_clean.pdb'
pdb2 		= parsePDB( inPDB ).select( 'protein' ).copy()
writePDB( inPDB_clean, pdb2 )

##


## Make span file from octopus
oSpan 		= os.path.join( oDir, inPDB[:-4] + '.span' )
oOctF 		= open( oSpan , 'w' )
cmdOCT		= [ 'perl', rosiOctPrl, octF ]

out = sp.call( cmdOCT, stdout=oOctF )
oOctF.close()
##

### Relax

cmd = [  rosiScrps, 
'-parser:protocol', protocolPth, 			# Path to Rosetta script (see above)
'-in:file:s', inPDB_clean,							# Input PDB structure
'-nstruct', '1', 							# Generate 1000 models
'-mp:setup:spanfiles', oSpan,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix', oDir,
'-out:overwrite',
'-packing:pack_missing_sidechains', 'false' ]

print 
print cmd
print 

sp.call( cmd )

###

## View in membrane
#rosiMemPyMol


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