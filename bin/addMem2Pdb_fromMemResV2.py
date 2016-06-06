## Marco Mravic DeGrado Lab April 2016
## 

# input 1: Input file,  
# input 2: Outfile file, pdb with DUM atoms as layers

##  python ~/CHAMP/bin/addMem2Pdb_fromMemRes.py match_1_0505.pdb > match_1_0505_MEM.pdb

import sys, os, numpy as np
from itertools import product
# parse input

txt = ''


help = '''
HETATM 1104 THKN MEM Y  64      18.235  -0.341   0.267  1.00  0.00           X  
HETATM 1105 CNTR MEM Y  64       3.295   0.968   0.005  1.00  0.00           X  
HETATM 1106 NORM MEM Y  64       3.255   0.708   0.970  1.00  0.00           X  

HETATM 4093  N   DUM  4093      22.000   6.000 -17.400                          
HETATM 4093  O   DUM  4093      22.000   6.000  17.400 
'''

memDict = {}
coords  = []
with open( sys.argv[1] ) as fin:
	for i in fin:

		if i[:4] 	=='ATOM':
			txt += i
			coords.append( np.array( [ float( i[30:38] ), float( i[38:46] ), float( i[46:54] ) ] ) )
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

#print unitNorm, thickness, np.linalg.norm( unitNorm ), d

point_out 		= unitNorm * 15
normal_out		= point_out + unitNorm
d_out = np.dot( -1* point_out, normal_out )

point_in 		= unitNorm * -15
normal_in		= point_in + unitNorm
d_in = np.dot( -1* point_in, normal_in )

#print d_out, normal_out

span = np.arange( -25.0, 25.0, 2 ) 

lines = ''

step = 4000
for k in product( span , repeat = 2 )  :
	
	# given x and y, solve for z with plane equation
	x, y = k[0], k[1]

	z_o =   round( ( -1 * d_out - normal_out[0] * x - normal_out[1] * y ) / normal_out[2], 3 )
	z_i =   round( ( -1 * d_in - normal_in[0] * x - normal_in[1] * y ) / normal_in[2], 3 )


	x += center[0]
	y += center[1]
	z_o += center[2]
	z_i += center[2]


	# ensure none are overlapping with protein coords
	skip_o = 0
	skip_i = 0
	for c in coords:

		if np.linalg.norm(  c - [ x, y, z_o ]  ) < 3:
			skip_o = 1
			break 

		if np.linalg.norm(  c - [ x, y, z_i ]  ) < 3:
			skip_i = 1
			break 
 

	if not skip_o:
		txt += 'HETATM %4d  O   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_o )
		step +=1

	if not skip_i:
		txt += 'HETATM %4d  N   DUM  %4d    %8.3f%8.3f%8.3f \n' % ( step, step, x,y,z_i )
		step +=1 

print txt
#print help
