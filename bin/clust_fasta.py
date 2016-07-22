# clean up redundant fasta files
# cluster sequences with appropriate BLOSOM/PAM matrix, average linkage
# python ~/CHAMP/bin/clust_fasta.py ~/CHAMP/EGFR/bbDesign_CTERMrnd2/match_1/redund_m1.fasta ~/CHAMP/EGFR/bbDesign_CTERMrnd2/match_1/redund_m1-dMat.pkl

import sys, os, numpy as np 

class Sequence:

	def __init__(self,seq, model, ps, score):
		self.seq 	= seq
		self.model 	= model
		self.ps 	= float( ps )
		self.score	= float( score ) 
		self.arrayIn= 10000
		self.clusLab= 10000

	def __repr__(self):
		return '#%s %.3f %.1f#' % ( self.model, self.ps, self.score )

unique = {}


## Read in file
with open( sys.argv[1] ) as fin:
	
	for i in fin:

		if i[0] == '>':
			model, score, ps = tuple( i.split()[1:] )
			continue
		if len(i.rstrip()) <1: continue
		if len(i) > 20: 
			seq = i.rstrip()

			try:
				if unique[seq].ps < ps:
					unique[seq] = Sequence( seq, model, ps, score )
				else: continue
			except KeyError:
				unique[seq] = Sequence(seq, model, ps, score)

## add index in matrix to Sequence object for lcustering aid later
stp = 0
for k in unique.keys():
	unique[k].arrayIn = stp
	stp += 1



from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
# choose matrix from a dry alignment run to calculate mean difference withn cluster
dmatrix = matlist.blosum85

#for a in pairwise2.align.globaldx("KEVLA", "EVL", dmatrix):

def fillMatrix( unique ):

	print '\nwriting distance matrix\n'
	matSize = len( unique.keys() )
	dMat 	= np.zeros( (matSize , matSize )  ) 
# calc all by all distance matrix.... first find avrg % seq similarity 
	r = 0
	c = 0

	for i in unique.keys():
		c = 0
		for j in unique.keys():

		#if i == j: continue
		#print pairwise2.align.globalds(i, j, dmatrix,-10, -10)[0]
			result = pairwise2.align.globalds(i, j, dmatrix,-10, -10)[0][2]
			dMat[r][c] = result
			c += 1
		r += 1
	return dMat/np.max( dMat )


#### Write the distance matrix for clustering to pickle file once only
import cPickle as pic
pklPath = sys.argv[2]
if not os.path.exists( pklPath ):
	dMat = fillMatrix( unique )
	pic.dump( dMat, open( pklPath, 'wb' ) )
else:
	dMat = pic.load( open( pklPath, 'rb' ) )


###### CLUSTERING #####

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, cophenet
import time
ts = time.time()
from collections import defaultdict

linkage_mat = linkage( dMat, method='complete')
for k in np.arange(1.0, 1.2, 0.002):

	clust 		= fcluster( linkage_mat, k )

	#print np.max(clust), k
	#print clust
	# organize clusters
	clusters 	= defaultdict(list)
	stp = 0
	for s in unique.keys():

		clusters[ clust[stp] ].append( unique[s] )
		unique[s].clusLab = clust[stp]

		stp +=1
	for k, v in clusters.items():
		repSeq = sorted( v, key= lambda x: x.ps , reverse=True)[0]
		print k, 'n=%d' % (len(v)), repSeq.ps, repSeq.seq

	print



	sys.exit()

#dendrogram( linkage_mat )












#import operator

#for v in sorted( unique.values(), key=operator.attrgetter('ps'), reverse=True):
#	print v, v.seq



