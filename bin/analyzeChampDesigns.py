# Marco Mravic DeGrado lab ucsf 2016 April

import sys, os, numpy as np
from prody import *
from PDButil import UnNatAA


##
class Model:
		def __init__(name, seq, score, uhb):
			self.name		= name
			self.seq		= seq
			self.score		= score
			self.uhb		= uhb


##

# input directory full of models
# parse for sequence of chain X (CHAMP), rosetta score, packing statistics, and buried hydrogen bonding moeties unpaired
# can also hard code this model's heptad assignment to note interfacial positions

# commadnline example
# python ~/CHAMP/bin/analyzeChampDesigns.py ~/CHAMP/a5/bb2Design/match_4/outputs/ > a5_match_4.fasta


# model starting heptad
hep = ''


models, seqs, scores, uhbs, packstats = [], [], [], [], []


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

			if i[:2] == 'uh':
				uhb = float( i.split()[1] )

	models.append( label )
	seqs.append( cSeq )
	scores.append( score )
	uhbs.append( uhb )
	packstats.append( packstat )

	# only read out subset 
	#if packstat < 0.50 or score > -210: continue

	print '>', label, score, packstat
	print cSeq, '\n'
#	print label, cSeq, score, packstat, uhb


packstats 	= np.array( packstats )
scores		= np.array( scores )


# plotting stuff
import matplotlib.pyplot as plt
#print len(cSeq)


plt.scatter( scores, packstats, color='blue',s=5,edgecolor='none')
plt.xlabel( 'Rosetta Score (REU)' )
plt.ylabel( 'Packing/Holes Score ' )
plt.title( sys.argv[2]  )
#plt.set_aspect(1./ax1.get_data_ratio())
plt.show()