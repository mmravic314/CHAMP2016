# Marco Mravic DeGrado Lab UCSF Biophysics Nov 2015
## hard coded search for alpha 5 integran matches in the 

# input 1: Directory of TM pairs, original
# input 2: ouput directory path
# input 3: tagrte helix to align the template helix to, given known PARING OF INPUT SEQUENCE TO TARGET HELIX STRUCTURE
# input 4: directory of EXTENDED TM pairs

## Example command line
# python ~/CHAMP/bin/templatePairsScreen_a5.py ~/CHAMP/Cluster-004/ ~/CHAMP/a5/bbHits/ ~/CHAMP/a5/RosMem_Target-input/a5_OPMout_clean_0001.pdb ~/CHAMP/Cluster-004_ext/ EGSYGVPLWIIILAILFGLLLLGLLIYILYKLG 13 20 [ASTG]\w\w[FMYWLIKR][ASG]\w\w[ILVFNQDEH]





import sys, os, re, numpy as np
from collections import defaultdict
from prody import *
from PDButil import *

# directory to deposit hits to
oDir = sys.argv[2]
if not os.path.exists(oDir):
	os.mkdir( oDir )


# Target sequence and corresponding indexing to string
#targ = 'EGSYGVPLWIIILAILFGLLLLGLLIYILYKLG'
#motifSrt, motifEnd =  13,20
targ = sys.argv[5]
motifSrt, motifEnd = int( sys.argv[6] ), int( sys.argv[7] )
motif_regex = sys.argv[8] 


## The membrane oriented helix to align the backbones of the sequence to
targHelix = parsePDB( sys.argv[3] )
targAtoms = targHelix.select( 'bb resnum %d to %d' % ( motifSrt + 1, motifEnd +1 ) ).copy()
#print [ x.getResname() for x in targAtoms.iterResidues() ]


# Save prody dimers after transformation to this found dict, titled by naming
foundDict=[]

for f in sorted( os.listdir( sys.argv[1] ) ):
	if f[-4:] != '.pdb': continue
	flg = 0

	print "Entering", f
	# seq seqs of each chain
	inPdb 	= parsePDB( os.path.join( sys.argv[1], f ), subset = 'bb' )
	seqA 	= ''.join( [ UnNatAA[x.getResname()] for x in inPdb.select( 'chain A' ).copy().iterResidues() ] )
	seqB	= ''.join( [ UnNatAA[x.getResname()] for x in inPdb.select( 'chain B' ).copy().iterResidues() ] )
	seqHash = { 'A':seqA, 'B':seqB }
	

	extPdb 	= parsePDB( os.path.join( sys.argv[4], f[:-4] + '_ext.pdb' ), subset = 'bb' )
	seqAext = ''.join( [ UnNatAA[x.getResname()] for x in extPdb.select( 'chain A' ).copy().iterResidues() ] )
	seqBext	= ''.join( [ UnNatAA[x.getResname()] for x in extPdb.select( 'chain B' ).copy().iterResidues() ] )
	extHash = { 'A':seqAext, 'B':seqBext }

	#print seqA
	#print seqB

	# search seqs of each chain for acceptable motifs by regex
	for k,v in seqHash.items():

		idN = f[:-4] + k

#		match = re.search(r'[ASG][AST]\w\w[ASTG][ASG]\w\w[ILV]')			## Alpha-2
		
		# find matching strings 
		matches = re.findall(r'%s' % motif_regex, v)			## alpha-5
		print matches, v

		if not matches: continue


		for match in matches:
			
			# locate this matching string to that on the extended pair
			match2 		= re.search( r'%s' % match, extHash[k] )
			srt, end 	= match2.start(), match2.end()
			print 'MATCH!', idN, match2.group(), srt, end, extHash[k]
		#	print extHash[k][srt:end], 
		#	print


			# move this pair of helices to align to membrane oriented helix (Input 3) at targeted motif site
			mobile 		= extPdb.select( 'bb chain %s resnum %d to %d' % (k, srt + 1, end  ) ).copy()
		#	print [ x.getResname() for x in mobile.iterResidues() ]

			rmsdB4 		= calcRMSD(mobile, targAtoms)
			rmsd, trMat = superpose( mobile, targAtoms )
			outPdb		= applyTransformation( trMat, extPdb.copy() )


			# Do rough distance cut off between opposing helix and template-matching helix... check for potential interaction or not
			ch2 	= 'B' 
			if k != 'A':
				ch2 = 'A'

			tempMtxhH 	= outPdb.select( 'bb chain %s resnum %d to %d' % (k, srt + 1, end  ) ).copy()
			oppoH 		= outPdb.select( 'bb chain %s ' % ch2 ).copy()

			dMat 		= buildDistMatrix( tempMtxhH, oppoH )
			# Skip if closest approach isn't in line with this motif interaction geometry: 4.3 cut off
			if dMat.min() > 4.3: 
				continue

			## Re-number chains: X for opposing/CHAMP helix, A for template-matching Helix
			chSwap = { k:'A', ch2:'X' }

			for i in outPdb.iterResidues():
				i.setChids(  chSwap[ i.getChid() ] )


			### Section Here need to be added for trimming the extended helical portions to TM only
#
#
#
			###



			## Save the residue range of the target motif in the Title of PDB
			title = idN + ' chain A %d to %d' % ( srt + 1, end )
			outPdb.setTitle( title )
			foundDict.append( outPdb )
			print title

			sys.exit()
			
	print


print 'found', len( foundDict ), 'eligible matches... now printied re-oriented PDBs to file'

# Print all these stored atom groups to file
mNum = 1
for i in foundDict:
	
	oPath = os.path.join( oDir, 'match_%d.pdb' % mNum )
	writePDB( oPath, i )
	mNum += 1
