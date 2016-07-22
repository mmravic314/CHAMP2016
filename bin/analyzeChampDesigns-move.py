# Marco Mravic DeGrado lab ucsf 2016 April

import sys, os, numpy as np
from prody import *
from PDButil import UnNatAA


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

# commadnline example
# python ~/CHAMP/bin/analyzeChampDesigns-move.py /home/xray/CHAMP/mouseEpoR/serZipper_30bb/serZip_30AA/outputs/ /home/xray/CHAMP/mouseEpoR/serZipper_30bb/serZip_30AA/Model_Validation/

# model starting heptad
hep = ''


models, seqs, scores, uhbs, packstats = [], [], [], [], []

Designs = {}

for m in os.listdir( sys.argv[1] ):
	if m[-4:] != '.pdb': continue

	mPath = os.path.join( sys.argv[1], m )

	label 	= os.path.splitext( m )[0].split('_')[-1]

	champ 	= parsePDB( mPath, chain='X' )
	cSeq 	= ''.join(  [ UnNatAA[ x.getResname() ] for x in champ.select('ca').iterAtoms() ] )



	# parse for energy, packstat, uhb
	with open( mPath ) as fin:
		for i in fin:
			if i[:2] == 'pa':
				packstat = float( i.split()[1] )

			if i[:2] == 'sc':
				score = float( i.split()[1] )

			#if i[:2] == 'uh':
			#	uhb = float( i.split()[1] )

#	models.append( label )
#	seqs.append( cSeq )
#	scores.append( score )
#	uhbs.append( uhb )
#	packstats.append( packstat )

	try:
		if Designs[ cSeq ].ps < packstat:
			Designs[ cSeq ] = Model( label, cSeq, score, packstat, mPath )
	except KeyError:
		Designs[ cSeq ] = Model( label, cSeq, score, packstat, mPath )

	# only read out subset 
	#if packstat < 0.50 or score > -210: continue

#	print '>', label, score, packstat
#	print cSeq, '\n'
#	print label, cSeq, score, packstat, uhb


#packstats 	= np.array( packstats )
#scores		= np.array( scores )

designs = []

## Top  30 packed
for k,v in sorted( Designs.items(), key=lambda x: x[1].ps, reverse=True )[:30]:
	designs.append( v.filePath )

# of remaining, top 30 high scoring
for k,v in sorted( Designs.items(), key=lambda x: x[1].score, reverse=True )[:60]:
	if len( designs ) == 60: break

	if v.filePath not in designs:
		designs.append( v.filePath )

# of remaining, 30 "randomly"
for k,v in Designs.items()[:90]:
	if len( designs ) == 90: break

	if v.filePath not in designs:
		designs.append( v.filePath )

targDir = sys.argv[2]
if not os.path.exists( targDir ):
	os.mkdir( targDir )

import subprocess as sp

for f in designs:
	print f
	sp.call( [ 'cp', f, targDir ] )


print 'Unique Designs:', len( Designs.keys() )


sys.exit()


# plotting stuff
import matplotlib.pyplot as plt
#print len(cSeq)


plt.scatter( scores, packstats, color='blue',s=5,edgecolor='none')
plt.xlabel( 'Rosetta Score (REU)' )
plt.ylabel( 'Packing/Holes Score ' )
plt.title( sys.argv[2]  )
#plt.set_aspect(1./ax1.get_data_ratio())
plt.show()