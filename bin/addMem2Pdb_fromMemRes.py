## Marco Mravic DeGrado Lab April 2016
## 

# input 1: Input file,  
# input 2: Outfile file, pdb with DUM atoms as layers

import sys, os, numpy as np
from itertools import product
# parse input

txt = ''

'''
HETATM 1104 THKN MEM Y  64      18.235  -0.341   0.267  1.00  0.00           X  
HETATM 1105 CNTR MEM Y  64       3.295   0.968   0.005  1.00  0.00           X  
HETATM 1106 NORM MEM Y  64       3.255   0.708   0.970  1.00  0.00           X  

HETATM 4093  N   DUM  4093      22.000   6.000 -17.400                          
HETATM 4093  O   DUM  4093      22.000   6.000  17.400 
'''

memDict = {}
with open( sys.argv[1] ) as fin:
	for i in fin:

		if i[:4] 	=='ATOM':
			txt += i
		elif i[:6] 	== 'REMARK':
			txt += i 
#		elif i[:3] 	== 'TER':
#			txt += i
		elif i[:6]  == 'HETATM':
			memDict[ i.split()[2] ] = np.array( [ float( i[30:38] ), float( i[38:46] ), float( i[46:54] ) ] )
		else:
			pass

#print txt

# calculate unit normal, 

center 		= memDict[ 'CNTR' ]
unitNorm 	= memDict[ 'NORM' ]	- center
thickness	= memDict[ 'THKN' ] - center
d = np.dot( unitNorm, thickness )

print unitNorm, thickness, np.linalg.norm( unitNorm ), d

d = 

span = np.arange( -20.0, 21.0, 2 ) 


for k in product( span, repeat=2 ):
	inner = np.array( [ k[0], k[1],  15.0 ] )
	outer = np.array( [ k[0], k[1], -15.0 ] )

	inner += center
