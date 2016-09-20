# clean up redundant fasta files
# cluster sequences with appropriate BLOSOM/PAM matrix, average linkage
# python ~/CHAMP/bin/clust_fasta.py ~/CHAMP/EGFR/bbDesign_CTERMrnd2/match_1/redund_m1.fasta ~/CHAMP/EGFR/bbDesign_CTERMrnd2/match_1/redund_m1-dMat.pkl

# 

import sys, os, numpy as np 

from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


#################### Global definitions, classes, functions  ####

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

dmatrix_HASH  = {}
dmatrix_HASH[30] = matlist.blosum30; dmatrix_HASH[35] = matlist.blosum35; dmatrix_HASH[40] = matlist.blosum40
dmatrix_HASH[45] = matlist.blosum45; dmatrix_HASH[50] = matlist.blosum50; dmatrix_HASH[55] = matlist.blosum55
dmatrix_HASH[60] = matlist.blosum60; dmatrix_HASH[65] = matlist.blosum65; dmatrix_HASH[70] = matlist.blosum70
dmatrix_HASH[75] = matlist.blosum75; dmatrix_HASH[80] = matlist.blosum80; dmatrix_HASH[85] = matlist.blosum85
dmatrix_HASH[90] = matlist.blosum90; dmatrix_HASH[95] = matlist.blosum95; dmatrix_HASH[100] = matlist.blosum100



#################### Global defs end  ###########################


#################  MAIN  ########################################

## Read in file
unique = {}
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



# choose matrix from a dry alignment run to calculate mean difference withn cluster (to the nearest 5)

sys.exit()


dmatrix = matlist.blosum85

#for a in pairwise2.align.globaldx("KEVLA", "EVL", dmatrix):

def fillMatrix( unique, dMat_key ):

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




