#calculate and plot association energy for a sequence

import sys, os, numpy as np
from prody import *
from PDButil import UnNatAA


# xray@dgl-xray:~/CHAMP/a5/m1_association$ python ~/CHAMP/bin/association_analysis.py  ~/CHAMP/a5/m1_association/ > ../a_scores.txt


## I/O

topDir = sys.argv[1]

##

model = []
score = []
packs = []

#header
print 'Model_ID delta_Rosetta_Energy Packstat_Score Sequence'

for model in os.listdir( topDir ):

	flg = 0



	if 'pdb' in model: continue
	midDir = os.path.join( topDir, model )
	sep_dir= os.path.join( midDir, 's_outputs/' )
	tog_dir= os.path.join( midDir, 't_outputs/' )

	# find scores of separated models
	s_score = []
	for pdb in os.listdir( sep_dir ):

		if not flg:

			inPDB 	= parsePDB( os.path.join( sep_dir, pdb ), subset='ca', chain='X' )
			seq 	= ''.join( UnNatAA[x] for x in inPDB.getResnames() )
			flg += 1


		with open( os.path.join( sep_dir, pdb ) ) as fin:
			for i in fin:

				if i[:6] == 'score_':
					s_score.append( float( i.split()[1] ) )

	# find scores of associated models 
	t_score = []
	ps_arr 	= []
	for pdb in os.listdir( tog_dir ):
		with open( os.path.join( tog_dir, pdb ) ) as fin:
			for i in fin:

				if i[:6] == 'score_':
					t_score.append( float( i.split()[1] ) )

				if i[:4] == 'pack':
					
					ps_arr.append( float( i.split()[1] ) )

	ps = max( ps_arr )
	a_sc = min( t_score ) - min( s_score )
	score.append( a_sc )
	packs.append(ps)
	print model, '%0.1f' % a_sc, '%0.3f' % ps, seq

packstats 	= np.array( packs )
scores		= np.array( score )


# plotting stuff
import matplotlib.pyplot as plt
#print len(cSeq)


plt.scatter( scores, packstats, color='blue',s=5,edgecolor='none')
plt.xlabel( 'Delta Rosetta Score (dREU)' )
plt.ylabel( 'Packing/Holes Score ' )
plt.title( 'Anti-Integin a5 CHAMP Association Energy, Unique Sequences (n=159)'  )
#plt.set_aspect(1./ax1.get_data_ratio())
plt.show()

