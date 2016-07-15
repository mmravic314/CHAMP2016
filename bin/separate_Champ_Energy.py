# Marco Mravic DeGrado lab ucsf 2016 July

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

# input directory full of models
# parse for sequence of chain X (CHAMP), rosetta score, packing statistics, and buried hydrogen bonding moeties unpaired
# can also hard code this model's heptad assignment to note interfacial positions

# Calculate association energy score... this should be done by re-doing design with helices apart, helices together... both with non-interfacial positions fixed
# This should set up 3 simulations... 

# commandline example
# python ~/CHAMP/bin/separate_Champ_Energy.py ~/CHAMP/a5/bb2Design/match_1/outputs/  ~/CHAMP/a5/m1_association/m1_separate.pdb  ~/rosetta/ ~/CHAMP/bin/helix_Relax-2Body.xml

# hard coded interfacial position index for non-interfacial residues... also ignore the last 4
# -- LLM LLL LLL AAL LGL LLG LLL VLL WKA NTR NK
# -- XXM XXX LXX AAX XGL XLG LXL VXX WXX XXX XX
# KK LLM LFL LFL AAI LGL ILG LFL VFL WKK
ignore 	= [ 0, 1, 3, 4 ,5, 7,8,11,12,15,19,22,23 ] # 25,26 LYS
# cut final residues off, exceppeting #28
rnd_AA 	= ['A', 'V', 'I', 'F', 'L', 'L','L','L','L','L']
newAA 	= {}
for x in ignore:
	newAA[x] = rnd_AA[ randint(0,9) ]

#[ rnd_AA[ randint(0,9) ] for x in ignore ]

h_len = 27
#models, seqs, scores, uhbs, packstats = [], [], [], [], []
Designs = {}

newFilesDir = os.path.dirname( sys.argv[2] )
spreadPDB 	= sys.argv[2]
rosiBase 	= sys.argv[3]
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')

repackXML 	= sys.argv[4]

timePer = []


### FIRST PASS
## for each design, parse output file for packstat & score... 
# remove deuplicate sequences by taking only best design (highest pack stat)
print '\nRecording Designs'
for m in os.listdir( sys.argv[1] ):
	if m[-4:] != '.pdb': continue

	mPath = os.path.join( sys.argv[1], m )
	label 	= os.path.splitext( m )[0].split('_')[-1]

	struct 	=  parsePDB( mPath )
	champ   = struct.select( 'chain X' ).copy()
	cSeq 	= ''.join(  [ UnNatAA[ x.getResname() ] for x in champ.select('ca').iterAtoms() ] )

	newSeq		 = ''
	interfaceSeq = ''
	for a in np.arange(25):
		if a in ignore:
			newSeq 		 += newAA[a]
		else:
			newSeq 		 += cSeq[a]
			interfaceSeq += cSeq[a]
		

	newSeq += 'KK'
	

	with open( mPath ) as fin:
		for i in fin:
			if i[:2] == 'pa':
				packstat = float( i.split()[1] )

			if i[:2] == 'sc':
				score = float( i.split()[1] )
	try:
		if Designs[ newSeq ].ps < packstat:
			Designs[ newSeq ] = Model( label, newSeq, score, packstat, mPath )
	except KeyError:
		Designs[ newSeq ] = Model( label, newSeq, score, packstat, mPath )

print 'Unique Designs Detected:', len( Designs.keys() )
# This double paring is super inefficient... but i dont want to rewrite this since Design is the rate limiting step
###

### SECOND PASS
## Go through designs again now... but skip those modesl whose labels arent the representative for that sequence


for m in os.listdir( sys.argv[1] ):
	if m[-4:] != '.pdb': continue

	start = time.time()

	### Prepare files for rosetta runs

	mPath = os.path.join( sys.argv[1], m )
	label 	= os.path.splitext( m )[0].split('_')[-1]

	struct 	=  parsePDB( mPath )
	champ   = struct.select( 'chain X' ).copy()
	cSeq 	= ''.join(  [ UnNatAA[ x.getResname() ] for x in champ.select('ca').iterAtoms() ] )

	newSeq		 = ''
	interfaceSeq = ''
	for a in np.arange(25):
		if a in ignore:
			newSeq 		 += newAA[a]
		else:
			newSeq 		 += cSeq[a]
			interfaceSeq += cSeq[a]
		

	newSeq += 'KK'

	## for each design, parse output file for packstat & score... 
	# remove deuplicate sequences by taking only best design (highest pack stat)

	model = Designs[newSeq]
	if model.name != label: 
		print 'model %s duplicate of %s... skip' % (label, model.name)
		continue

	##


	toRepack 	= struct.select( 'resnum 1 to 58' ).copy()
	innerDir 	= os.path.join( newFilesDir , os.path.basename( mPath[:-4] ) + '/' )
	print 'entering...', innerDir, newSeq

	if not os.path.exists(innerDir): 

		os.mkdir( innerDir )
		

	resF 		= os.path.join( innerDir, label + '.resfile' )
	Xcnt 		= 0
	txt = 'START\n'
	for r in toRepack.iterResidues():

		if r.getChid() == 'A':
			#print r.getResname(), r.getResnum(), 'A'
			txt += '%d %s PIKAA %s\n' % ( r.getResnum(), 'A', UnNatAA[ r.getResname() ] )

		else:
			#print newSeq[Xcnt], r.getResnum(), 'X'
			txt += '%d %s PIKAA %s\n' % ( r.getResnum(), 'X', newSeq[Xcnt] )
			Xcnt += 1
	outRes = open( resF, 'w' )
	outRes.write( txt )
	outRes.close()


	newPath 	= os.path.join( innerDir, os.path.basename( mPath) )
	writePDB( newPath, toRepack )

	### End rosetta file prep

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
		sp.call( cmdSpan, stdout=FNULL, stderr=sp.STDOUT  )

	# run separates 10x
	sOutDir = os.path.join( innerDir, 's_outputs/' )
	n = '10'
	if not os.path.exists( sOutDir ):
		os.mkdir( sOutDir )

	s_cmd 	= [ rosiScrps, 
'-parser:protocol', repackXML, 					# Path to Rosetta script (see above)
'-in:file:s', spreadPDB,						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:setup:spanfiles', spanPath,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix',sOutDir,
'-packing:resfile', resF,
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

	sp.call(  s_cmd , stdout=FNULL, stderr=sp.STDOUT  )


	# run together  10x 
	tOutDir = os.path.join( innerDir, 't_outputs/' )
	n = '10'
	if not os.path.exists( tOutDir ):
		os.mkdir( tOutDir )



	t_cmd 	= [ rosiScrps, 
'-parser:protocol', repackXML, 					# Path to Rosetta script (see above)
'-in:file:s', newPath,						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:setup:spanfiles', spanPath,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix',tOutDir,
'-packing:resfile', resF,
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

	sp.call(  t_cmd, stdout=FNULL, stderr=sp.STDOUT  )



	############
	elapsed = time.time() - start
	timePer.append( elapsed )
	print 'This took', elapsed, 'sec', '\n'


	#sys.exit()
	continue

	## Given the sequence for this model, put it into a new file, thread this sequence and calculate
	## sequence is 32 to 58


	# parse for energy, packstat, uhb
	with open( mPath ) as fin:
		for i in fin:
			if i[:2] == 'pa':
				packstat = float( i.split()[1] )

			if i[:2] == 'sc':
				score = float( i.split()[1] )


	try:
		if Designs[ newSeq ].ps < packstat:
			Designs[ newSeq ] = Model( label, cSeq, score, packstat, mPath )
	except KeyError:
		Designs[ newSeq ] = Model( label, cSeq, score, packstat, mPath )


#print 'Unique Designs:', len( Designs.keys() )



sys.exit()

