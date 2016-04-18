# Marco Mravic DeGrado Lab UCSF Biophysics Nov 2015
## hard coded search for alpha 5 integran matches in the 


targ = 'PLWIIILAILFGLLLLGLLIYVLYK'
motif= ''

import sys, os, re
from collections import defaultdict
from prody import *
from PDButil import *


foundDict={}

for f in sorted( os.listdir( sys.argv[1] ) ):
	if f[-4:] != '.pdb': continue
	flg = 0

	print "Entering", f
	# seq seqs of each chain
	inPdb 	= parsePDB( os.path.join( sys.argv[1], f ), subset = 'ca' )
	seqA 	= ''.join( [ natAA[x] for x in inPdb.select( 'chain A' ).getResnames() ] )
	seqB	= ''.join( [ natAA[x] for x in inPdb.select( 'chain B' ).getResnames() ] )
	seqHash = { 'A':seqA, 'B':seqB }

	#print seqA
	#print seqB

	# search seqs of each chain for acceptable motifs by regex
	for k,v in seqHash.items():

		idN = f[:-4] + k

#		match = re.search(r'[ASG][AST]\w\w[ASTG][ASG]\w\w[ILV]')			## Alpha-2
		match = re.search(r'[ASTG]\w\w[FMYWLIKR][ASG]\w\w[ILVFNQDEH]', v)			## alpha-5
		if match:
			print 'MATCH!', idN, v, match.group()
			#sys.exit()

	print

	






















sys.exit()



# Check input arguments
if len(sys.argv) != 5:
	print "ERROR: not enough input arguments specified"
	print	help
	sys.exit()


# Parse input sequence file
tar_seq	= ''		# target (biological) sequence
tar_mot = ''		# interaction motif to utilize
space 	= 0
spacer	= ''
spaceRE	= ''
tar_re	= ''

with open( sys.argv[2] ) as file:
	for i in file:
		tar_seq += i.rstrip()[0]
		

		if i.rstrip()[-1] == '*':
			if len(tar_mot) == 0:
				tar_mot += i.rstrip()[0]
				tar_re  += tar_mot
				continue
			else:
				tar_mot += spacer  + i.rstrip()[0]
				tar_re  += spaceRE + i.rstrip()[0]
				continue

		if len(tar_mot) > 0:
			spacer  += 'X'
			spaceRE += '\w'


chainRv = { 'A':'B', 'B':'A' }			# return opposite chain of chain given

## Given a motif to match, namely Small-XXX-G-XXX-Large
#  check a few quantitative criteria for template pair PDB
#  Input: ProDy atom group, chain of match, tuple of indices motif aligned to, sequence of match)

def evalMatch( pdbObj, chID, spnTup, seq ):
	print "Evaluating potential match:", chID, seq, spnTup, '...',
	
	# Check if proline is within sequence -2 to +8 of match to ATS_XXX_G; skip if too short
	try:
			if 'P' in seq[ spnTup[0] - 2: spnTup[-1] + 3 ]:	
					print "Proline found in target region"
					return False
			if len(seq) < ( spnTup[-1] + 5 ):
					print "Template sequence too short"
					return False
	except IndexError:
					print "Template sequence too short"
					return False

	# Check if these two Ca's are within striking distance of other helix's interface (7.3 A)
	distMat = buildDistMatrix(   pdbObj.select( 'chain %s resnum %d %d' % ( chID, spnTup[0], spnTup[1]) ).getCoords()   ,  inPdb.select( 'chain %s' % ( chainRv[chID] ) ).getCoords() )
	if min( distMat[0] ) > 7.5 or min( distMat[1] ) > 6.7:
		print "Target motif matched not at interface in template PDB"
		return False

	print "\nMATCH ACCEPTED"
	return True


print "\nInput target sequence", tar_seq, '\nWith motif to align on template pairs', tar_mot, tar_re
print '\nEntering template pairs directory...\n'

# hack 
tar_re = r'[AST]\w\w\wG'
n_smal = ['A', 'S', 'T']
n2smal = ['A', 'G']
n_lrg  = ['L', 'I', 'F']

# Enter each template pair PDB, and do a sequence alignemnt to motif if found, 

for f in sorted( os.listdir( sys.argv[1] ) ):
	if f[-4:] != '.pdb': continue
	flg = 0

	print "Entering", f

	inPdb 	= parsePDB( os.path.join( sys.argv[1], f ), subset = 'ca' )
	seqA 	= ''.join( [ natAA[x] for x in inPdb.select( 'chain A' ).getResnames() ] )
	seqB	= ''.join( [ natAA[x] for x in inPdb.select( 'chain B' ).getResnames() ] )
	seqHash = { 'A':seqA, 'B':seqB }
	
	fnd = defaultdict(list)
	# Check motifs in chain A sequence, save tuples of index where motif aligns
	i = 0
	while (i+8) < len( seqA ) - 1:
		if seqA[i] in n_smal and seqA[i+4] in n2smal and seqA[i+8] in n_lrg: 
			fnd['A'].append( ( i, i+4, i+8 ) )
		i += 1
	# Check motifs in chain B sequence, as above ""
	i = 0
	while ( i+8 ) < len( seqB ) - 1:
		if seqB[i] in n_smal and seqB[i+4] in n2smal and seqB[i+8] in n_lrg: 
			fnd['B'].append( ( i, i+4, i+8 ) )
		i += 1

	# For each potential alignment, evaluate match by few critera (discussed in function description)
	# return boolean TRUE if acceptable for sequence threading
	matches = defaultdict(list)
	for chain, tup in fnd.items():
		for match in tup:
			if evalMatch( inPdb, chain, match, seqHash[chain] ):
				matches[chain].append( match )
				



	for k,v in matches.items():
		print 'chain', k, seqHash[k]
		print v
		print 

	print 


## just SXXXG search does not give nay good templates 
#  4ain-030_038-0138_0149_C-0394_0423_C.pdb CHAIN B SXXXG Does not align w/ interface
#  4a01-020_022-0186_0216_B-0287_0311_B.pdb CHAIN B SXXXG bcdab S not a interfacial positions (but okay)
#  3gia-000_006-0010_0036_A-0220_0245_A.pdb CHAIN A SXXXG bcdef S nor G at interface
#  3qf4-006_008-0036_0065_B-0148_0175_B.pdb CHAIN A Good alignment but proline downstream. Chain B has good SXXXT at interface?
#  3b9w-023_025-0042_0065_C-0104_0124_C.pdb CHAIN A Not in membrane or at interface
#  2nq2-002_006-0099_0113_A-0236_0256_A.pdb CHAIN A Not inter interaction interface
#  2jln-000_006-0029_0040_A-0246_0274_A.pdb CHAIN B SXXXG Does not align w/ interface
#  2cfq-005_009-0220_0249_A-0347_0371_A.pdb CHAIN B At interface but bottom of template helix

## [ATS]XXXG template search

#### Searched [ATS]XXX[GA]XXX[LIF] 10/11/15
# 4gc0-006_012-0277_0291_A-0404_0423_A.pdb B RGKALAIAVAANWLANYFVS (5, 9, 13)  BAD
# 3k3f-001_003-0046_0066_A-0097_0121_A.pdb B PAMLGYVALNGAFTTIIMASLLNFL (7, 11, 15) GOOD
# 3ayf-001_009-0248_0274_A-0568_0586_A.pdb B EVWIALGAVFSALEVIPLT (7, 11, 15)  C-term proline overwritable?
# 2w2e-017_020-0168_0191_C-0242_0266_C.pdb A RTRGLFLEAFGTAILCLTVLMLAV (8, 12, 16) Good! needs extension? THREADS GREAT!!
# 2f2b-010_013-0143_0163_B-0222_0243_B.pdb A YWNAMLAEVVGTFLLMITIMG (6, 10, 14) Doesn't interface
# 2f2b-000_002-0004_0029_A-0100_0122_A.pdb A LTKRCIAEFIGTFILVFFGAGSAAVT (6, 10, 14) Doesn't interface