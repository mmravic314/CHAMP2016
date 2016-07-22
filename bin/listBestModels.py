# listBestModels.py
#  python ~/CHAMP/bin/listBestModels.py ~/CHAMP/a5/m1_association/ ~/CHAMP/a5/model_list.txt

# input 1: directory at root of tree. contaning subdir for each model, which has subdir t* and s* for separate and spread designs
# input 2: path of output list to write the path of teh best associated final pose for each model to

import os, sys

## I/O ##

treeHead 	= sys.argv[1]
outPath 	= sys.argv[2]

## ##  ##

txt = ''
for subdir in os.listdir( treeHead ):

	path = os.path.join( treeHead, subdir )
	if not os.path.isdir( path ): continue

	# find lowest scoring separated model, for baseline of association
#	sepDir 	= os.path.join( path, 's_outputs' )
	score = []
#	for f in [ x for x in os.listdir(sepDir) if x[-3:] == 'pdb']:

#		with open( os.path.join( sepDir, f ) ) as fin:
#			for i in fin:

#				if i[:2] == 'sc':
#					score.append( float( i.split()[1] ) )

#	print min(score_s) 

	# find lowest scoring associated model, and remember it's 
	sepDir 	= os.path.join( path, 't_outputs' )
	score_s = []
	for f in [ x for x in os.listdir(sepDir) if x[-3:] == 'pdb']:

		with open( os.path.join( sepDir, f ) ) as fin:
			for i in fin:

				if i[:2] == 'sc':
					score.append( (float( i.split()[1] ), os.path.join( sepDir, f )  ) )
	txt += sorted( score )[0][1] + '\n'



#	sys.exit()
oF = open( outPath, 'w' )
oF.write(txt )




