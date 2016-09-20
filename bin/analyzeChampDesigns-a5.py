# Marco Mravic DeGrado lab ucsf 2016 April

# example command line 
# 
# python ~/CHAMP/bin/analyzeChampDesigns-a5.py ~/CHAMP/a5/bb2Design/match_1/outputs/ ~/CHAMP/a5/bb2Design/match_1/dmatrix.pkl 'CHAMP anti-a5 Backbone 1' > ~/CHAMP/a5/bb2Design/match_1/a5_m1_fasta_summary.txt
# 
# python ~/CHAMP/bin/analyzeChampDesigns-a5.py ~/CHAMP/a5/bb2Design/match_4/outputs/ ~/CHAMP/a5/bb2Design/match_4/dmatrix.pkl 'CHAMP anti-a5 Backbone 4' > ~/CHAMP/a5/bb2Design/match_4/a5_m4_fasta_summary.txt


import sys, os, numpy as np
from prody import *
from PDButil import UnNatAA

from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


#################
class Model:
		def __init__(self,name, seq, score, ps, filepath):
			self.name		= name
			self.seq		= seq
			self.score		= score
			self.ps 		= ps
#			self.filePath	= filepath
			self.cluster 	= 1000

		def __repr__( self ):
			return self.seq

def roundFive(x, base=5):
	return int(base * round(float(x)/base))


dmatrix_HASH  = {}
dmatrix_HASH[30] = matlist.blosum30; dmatrix_HASH[35] = matlist.blosum35; dmatrix_HASH[40] = matlist.blosum40
dmatrix_HASH[45] = matlist.blosum45; dmatrix_HASH[50] = matlist.blosum50; dmatrix_HASH[55] = matlist.blosum55
dmatrix_HASH[60] = matlist.blosum60; dmatrix_HASH[65] = matlist.blosum65; dmatrix_HASH[70] = matlist.blosum70
dmatrix_HASH[75] = matlist.blosum75; dmatrix_HASH[80] = matlist.blosum80; dmatrix_HASH[85] = matlist.blosum85
dmatrix_HASH[90] = matlist.blosum90; dmatrix_HASH[95] = matlist.blosum95; dmatrix_HASH[100] = matlist.blosum100




def fillMatrix( unique ):

	matSize = len( unique )
	print '\nwriting distance matrix, n = %d\n' %  matSize

	dMat 	= np.zeros( (matSize , matSize )  ) 
	s_len 	= len( unique[0] )


	# calc a miniBatch (1%) of the all by all similarity matrix to first find average % seq similarity 
	stop = 0
	seqID = []
	r = 0
	while stop < 0.01 * matSize**2 :
		for i in unique:
			c = 0
			for j in unique:

				if r == c: continue 	# Skip diagonal of matrix

				same 	= 0
				for a,b in zip( i, j ):
					if a == b: same += 1
				sID 	= float( same ) / s_len
				seqID.append( sID * 100 )
				stop += 1

				c += 1
			r = 1

	# choose the blosom similarity matrix to use
	seqID = np.array( seqID )
	level = roundFive( np.mean( seqID ) )
	blosom  = dmatrix_HASH[ level ]

	# then fill out the sequence pairwise matrix
	r = 0
	for i in unique:
		c = 0
		for j in unique:

			result 			= pairwise2.align.globalds(i, j, blosom, -10, -10 )[0][2]
			dMat[r][c] 		= result

			c += 1
		r += 1

	return dMat


rnd_AA 	= ['A', 'V', 'I', 'F', 'L', 'L','L','L','L','L']
#for x in ignore:
	#newAA[x] = rnd_AA[ randint(0,9) ]

#################

# input directory full of models
# parse for sequence of chain X (CHAMP), rosetta score, packing statistics, and buried hydrogen bonding moeties unpaired

### In this file version, hard code the interfacial positions & their forces identities

#####  match 4, '_' is designable, 'X' is non-interfacial/fixed, other positions manually set
### clip 0, KKW [W is 2], (residue 5 - 30 ), KK
### the count starts where residue 1 is gone, and residue 2 is index 0. index 2 set to 'W'
# PLL LLL LLI ILI MLI IAL IVA IIL MLL AIL LA
# -KK X__ XX_ XXX __X __X X__ XX_ _XX __X XAK K  <- pseduo-coded format of these designs
#  KK XLL XXI XXX MLX IAX XVA XXL MXX AIX XAK K  <- what the given design would look like, X's replaced by randoms
#  KK LFA FVL ILA VMV FIA LLV ALL LML FAI LKK 	 <- actual design 
fixed 	= 'KKLFAFVLILAVMVFIALLVALLLMLFAILKK'
ignore 	= np.array( [ 2, 5, 6, 8, 9, 10, 13, 16, 17, 20, 21, 24, 25, 28, 29 ] ) -1



## match 1  ##########

# hard coded interfacial position index for non-interfacial residues... also ignore the last 4
# -- LLM LLL LLL AAL LGL LLG LLL VLL WKA NTR NK
# -- XXM XXX LXX AAX XGL XLG LXL VXX WXX XXX XX 
# KK VLM VVL LLL AAL LGI LLG LIM VLW VKK		 <- actual design, m1_779 

## Hash this in for match 1
#ignore 	= [ 0, 1, 3, 4 ,5, 7,8,11,12,15,19,22,23 ] # 25,26 LYS
#fixed 	= 'KKVLMVVLLLLAALLGILLGLIMVLWVKK'

# cut final residues off, exceppeting #28




######################


# commadnline example
# python ~/CHAMP/bin/analyzeChampDesigns.py ~/CHAMP/a5/bb2Design/match_4/outputs/ > a5_match_4.fasta

#### STEP 1) parse sequence & energy data & link through model object. Take highest packstat for identical seqs 

models, seqs, scores, uhbs, packstats = [], [], [], [], []
Designs = {}

for m in os.listdir( sys.argv[1] ):
	if m[-4:] != '.pdb': continue

	mPath = os.path.join( sys.argv[1], m )

	label 	= os.path.splitext( m )[0].split('_')[-1]

	
	## mod to accomodate weird filename differences from cluster design versus local design
#	label 	= os.path.splitext( m )[0].split('out')[-1].split('_')[0]
	##


	champ 	= parsePDB( mPath, chain='X' )
	cSeq 	= ''.join(  [ UnNatAA[ x.getResname() ] for x in champ.select('ca').iterAtoms() ] )

	# replace the cSeq non-interfacials and modify as described above
	step 	= 0

	newSeq 	= 'KK'			# this is hardcoded for match 4
	for k in cSeq[2:-2]:	# this is hardcoded for match 4

#	newSeq 	= 'KK'			# this is hardcoded for match 1
#	for k in cSeq[:25]:		# this is hardcoded for match 1

		if step in ignore:
			newSeq += 'C'
			step += 1
			continue

		newSeq += k


		step += 1

	newSeq += 'KK'				# this is hardcoded for match 4 & 1

	# parse for energy, packstat from rosetta output file
	with open( mPath ) as fin:
		for i in fin:
			if i[:2] == 'pa':
				packstat = float( i.split()[1] )

			if i[:2] == 'sc':
				score = float( i.split()[1] )


	# only read out subset 
#	if packstat < 0.60 or score > -220: continue

	# put into model object 
	try:
		if Designs[ newSeq ].ps < packstat:
			Designs[ newSeq ] = Model( label, newSeq, score, packstat, mPath )
	except KeyError:
		Designs[ newSeq ] = Model( label, newSeq, score, packstat, mPath )


###### STEP 2) cluster the unique sequences, and add labels to the design objects 
# Here it's hierachical with separation at the level of the seqID, given blosom scores (Blosom80 for seq ID ~80%)

import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage, fcluster

#### Write the similarity matrix for clustering to pickle file once only
import cPickle as pic
pklPath = sys.argv[2]
if not os.path.exists( pklPath ):
	sMat = fillMatrix( Designs.keys() )
	pic.dump( sMat, open( pklPath, 'wb' ) )
else:
	sMat = pic.load( open( pklPath, 'rb' ) )

# calculate the average pair-wise non-identical similaity from the matrix:
data = []
step = 0
for k in sMat:
	data.extend( k[step + 1:] )
	data.extend( k[:step] )
	step += 1

# normalize the matrix by the maximum non-idential similarity blosum score
sMat = sMat/ np.max( data )


# fix those identical positions at 1... this is and the previous move are both hairy, break communitive property
for k in np.arange( len( Designs.keys() ) ):
	sMat[k][k] = 1.0

# convert similarity to distance matrix, then convert redundant n*n square matrix form into a condensed nC2 array
dMat = 1 - sMat
distArray = ssd.squareform( dMat )

# Clustering cutoff at the average distance between all members
### Here you can hardcode in 
cutoff 	  =	 np.mean( distArray )  

print 'complete linkage clustering, cutoff =', cutoff
# perform complete linkage clustering & return clusters at described cutoff
linkMat 	= linkage( distArray, method='complete', metric='euclidean')
clust 		= fcluster( linkMat, cutoff, criterion='distance')

print 'Unique clusters found:', len( set(clust) ), '\n'		## This was 33 at 1.0 angstrom RMSD

# assign clusters and write to file
for k, cl in zip( Designs.keys(), clust ):
	Designs[k].cluster = cl

## print out designs organized by cluster to stdout... also collect non-redundant scores & packstat
## Also collect those suggested to follow-up, based on best in cluster & 'okay' in general
prvModel = ''
print '# fasta Design_ID Rosetta_score Rosetta_Holes_Packstat Cluster_ID\n'
scores 		= []
packstats 	= []
suggested   = []
clus_models = []


for v in sorted( Designs.values(), key=lambda x: (x.cluster, x.ps) ):
	
	## store score & ps
	scores.append( v.score )
	packstats.append( v.ps )



	if prvModel == '':
		print '________  cluster 1  ____________\n' 
		prvModel = v
		

	# at stat of a new cluster, check if the last of previous is good enough to suggest, and print cluster header
	elif prvModel.cluster != v.cluster :
		print '________  cluster %s  ____________\n' % v.cluster

	# review the cluster and pick the 'best'
	# no TM tryptophans, highest packstat in cluster, ps greater than 0.5
		if len( clus_models ) > 0:
			suggested.append( sorted( clus_models, key=lambda x: x.ps, reverse=True )[0] )
			clus_models = []



	info = '> design %s %0.3f %0.3f %s' % ( v.name, v.score, v.ps, v.cluster )
	print info, '\n', v.seq, '\n'
	
	prvModel = v
	if 'W' not in prvModel.seq[6:-5] and prvModel.ps > 0.5: 
		clus_models.append( v ) 


# Handle the final model (best of last cluster possibly)
if len( clus_models ) > 0:
	suggested.append( sorted( clus_models, key=lambda x: x.ps, reverse=True )[0] )




print '\n', 'So-called best, diverse designs (picked automatically) '
for v in suggested:
	info = '> design %s %0.3f %0.3f %s' % ( v.name, v.score, v.ps, v.cluster )
	print info, '\n', v.seq, '\n'


# plotting stuff
import matplotlib.pyplot as plt

scores 		= np.array( scores )
packstats 	= np.array( packstats )

plt.scatter( scores, packstats, color='blue',s=5,edgecolor='none')
plt.xlabel( 'Rosetta Score (REU)' )
plt.ylabel( 'Packing/Holes Score ' )
plt.title( sys.argv[3]  )
plt.show()