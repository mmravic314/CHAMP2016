# Marco Mravic DeGrado lab ucsf 2016 July, EDIT just to work with a5 match 4 backbones... some hard coded things
# big change is taking an input pathlist file to the models to work with, and picking the fixed sequence

##
# python ~/CHAMP/bin/separate_Champ_Energy-a5m4.py ~/CHAMP/a5/bb2Design/match_4/model_listS.txt  ~/CHAMP/a5/m4_association/m4_separate.pdb  ~/rosetta/ ~/CHAMP/bin/helix_Relax-2Body.xml

import sys, os, numpy as np, subprocess as sp, time
from prody import *
from PDButil import UnNatAA
from random import randint



##
class Model:
		def __init__(self,name, seq, score, ps, filepath):
			self.name		= name
			self.seq		= seq
			self.score		= score
#			self.uhb		= uhb
			self.ps 		= ps
			self.filePath	= filepath


##

# input list of paths to modelss 
# parse for sequence of chain X (CHAMP), and replace with a fixed/random set of apolar groups

# Calculate association energy score... this should be done by re-doing design with helices apart, helices together... both with non-interfacial positions fixed
# This should set up 3 simulations... 


### In this file version, hard code the interfacial positions & their forces identities

#####  match 4, '_' is designable, 'X' is non-interfacial/fixed, other positions manually set
### clip 0, KKW [W is 2], (residue 5 - 30 ), KK
### the count starts where residue 1 is gone, and residue 2 is index 0. index 2 set to 'W'
# PLL LLL LLI ILI MLI IAL IVA IIL MLL AIL LA
# -KK X__ XX_ XXX __X __X X__ XX_ _XX __X XAK K  <- pseduo-coded format of these designs
#  KK XLL XXI XXX MLX IAX XVA XXL MXX AIX XAK K  <- what the given design would look like, X's replaced by randoms
#  KK LFA FVL ILA VMV FIA LLV ALL LML FAI LKK 	 <- actual design 
fixed 	= 'KKLFAFVLILAVMVFIALLVALLLMLFAILKK'
ignore 	= [ 1, 4, 5, 7, 8, 9, 12, 15, 16, 19, 20, 23, 24, 27, 28 ]
instead	= [ fixed[ k +2 ] for k in ignore ]



### I/O ### 

newFilesDir = os.path.dirname( sys.argv[2] )
spreadPDB 	= sys.argv[2]							# Make this file if it does not exist
rosiBase 	= sys.argv[3]
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')

repackXML 	= sys.argv[4]

timePer = []

###########



###### MAIN ######

## loop through input path list file, working on each
#  Modify the sequence and structure: shorten, add LYS to the end, & change non-interfacial. Consistent w/ input seperate models 
#  Prepare the pakcing resfile
##

with open( sys.argv[1] ) as fin:
  for m in fin:
	start = time.time()

	### Prepare files for rosetta runs

	mPath 	= m.rstrip()
	label 	= os.path.splitext( m )[0].split('_')[-1] 		# This needs to be flagged for design filenames done on a Cluster

	struct 	=  parsePDB( mPath )
	champ   = struct.select( 'chain X' ).copy()
	cSeq 	= ''.join(  [ UnNatAA[ x.getResname() ] for x in champ.select('ca').iterAtoms() ] )

# replace the cSeq non-interfacials and modify as described above
	step 	= 0
	newSeq 	= 'KK'			# this is hardcoded for match 4
	for k in cSeq[2:-2]:	# this is hardcoded for match 4

#	newSeq 	= 'KK'			# this is hardcoded for match 1
#	for k in cSeq[:25]:		# this is hardcoded for match 1

		if step in ignore:
			#newSeq += instead[ignore.index(step)]
			newSeq += 'x'
			step += 1
			continue

		newSeq += k
		step += 1

	newSeq += 'KK'				# this is hardcoded for match 4 & 1

	#sys.exit()

	## This sections to choose what residues of bnoth helices to select... here all are used.
#	champResSet	= ' '.join([ str(z) for z in struct.select( 'ca chain X' ).getResnums()[:] ])    
#	toRepack 	= struct.select( '(chain A) or ( chain X and resnum %s )' % champResSet ).copy()

	toRepack	= struct.copy()
	innerDir 	= os.path.join( newFilesDir , os.path.basename( mPath[:-4] ) + '/' )
	##

	print 'entering...', innerDir, newSeq
	if not os.path.exists(innerDir): 
		os.mkdir( innerDir )
	
	# Write resfile, renumber both helices to start from 1... Changes only if the structure was truncated/extended
	resF 		= os.path.join( innerDir, os.path.basename( mPath )[:-4] + '.resfile' )
	Xcnt 		= 0
	resN 		= 1
	txt = 'START\n'
	for r in toRepack.iterResidues():

		r.setResnum( resN )

		if r.getChid() == 'A':
			#print r.getResname(), r.getResnum(), 'A'
			txt += '%d %s PIKAA %s\n' % ( r.getResnum(), 'A', UnNatAA[ r.getResname() ] )

		else:
			if r.getResname() == 'MEM': continue		# Skip rosetta membrane residue
			txt += '%d %s PIKAA %s\n' % ( r.getResnum(), 'X', newSeq[ Xcnt ] )
			Xcnt += 1
		resN += 1

	outRes = open( resF, 'w' )
	outRes.write( txt )
	outRes.close()

	newPath 	= os.path.join( innerDir, os.path.basename( mPath) )
	writePDB( newPath, toRepack )

	### End rosetta file prep


	## Make the reference 'separated PDB ' if it does not already exist
	if not os.path.exists( spreadPDB ):

		## Use membrane normal vector & center to calculate  some vector that lies in membrane parallel, to spread helices parallel this vector
		# correct for normal vector's position to center, and use plane point as origin 
		#  norm(a,b,c) *dot* ( vect(x,y,z) - cent(0,0,0) ) = 0; solve for (x,y,z) given x & y arbitrary
		memCent, memNorm = toRepack.select( 'name CNTR' ).getCoords()[0], toRepack.select( 'name NORM' ).getCoords()[0]
		memNorm =  memNorm - memCent 

		x, y = 1., 1.
		z	 =  -1 * ( memNorm[0] * x + memNorm[1] * y ) / memNorm[2]

		vect = np.array( [x,y,z] )
		vect = vect / np.linalg.norm( vect ) 

		helix1 = toRepack.select( 'chain A' )
		helix2 = toRepack.select( 'chain X' )
		moveAtoms( helix1, by=-30 * vect )
		moveAtoms( helix2, by= 30 * vect )
		spread = helix1 + helix2 + toRepack.select( 'resname MEM' ) 
		writePDB( spreadPDB, spread )

	# end of making separated PDB


	# dump rosetta std output... noisey
	FNULL = open(os.devnull, 'w')

	# run to get the span file
	os.chdir( innerDir )
	spanPath = newPath[:-4] + '.span'
	if not os.path.exists( spanPath ):

		cmdSpan = [ rosiSpanGen, 
		'-database', rosiDB, 
		'-in:file:s', newPath
		]	
#		sp.call( cmdSpan, stdout=FNULL, stderr=sp.STDOUT  )

	# run separates 10x
	sOutDir = os.path.join( innerDir, 's_outputs/' )
	
	n = '10'
	
	if not os.path.exists( sOutDir ):
		os.mkdir( sOutDir )

	s_cmd 	= [ rosiScrps, 
'-parser:protocol', repackXML, 					# Path to Rosetta script (see above)
'-in:file:s', spreadPDB,						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix',sOutDir,
'-packing:resfile', resF,
'-out:overwrite',
'-mp::setup::spans_from_structure', '1', 	# Take span file from the membrane residue 
'-packing:pack_missing_sidechains', '0' ]

	sp.call(  s_cmd , stdout=FNULL, stderr=sp.STDOUT  )
#	sp.call(  s_cmd  )


	# run together  10x 
	tOutDir = os.path.join( innerDir, 't_outputs/' )
	n = '10'
	if not os.path.exists( tOutDir ):
		os.mkdir( tOutDir )


	t_cmd 	= [ rosiScrps, 
'-parser:protocol', repackXML, 					# Path to Rosetta script (see above)
'-in:file:s', newPath,						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix', tOutDir,
'-packing:resfile', resF,
'-mp::setup::spans_from_structure', '1', 	# Take span file from the membrane residue 
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

	sp.call(  t_cmd, stdout=FNULL, stderr=sp.STDOUT  )
#	sp.call( t_cmd )
#	sys.exit()

	############
	elapsed = time.time() - start
	timePer.append( elapsed )
	print 'This took', elapsed, 'sec', '\n'

print 'total (s):', sum( timePer ), np.mean( timePer )  

### Usually takes 2-3 minutes per model for the 20 trials 
