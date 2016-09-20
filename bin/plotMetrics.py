import sys, os, numpy as np
import matplotlib.pyplot as plt


vacRMSD = {}

with open( sys.argv[2] ) as fin:
	for f in fin:
		name, rmsd = tuple( f.split()[:2] )
		vacRMSD[name] = float( rmsd )

rmsd = []
a_sc = []
p_st = []

with open( sys.argv[1] ) as fin:
	for f in fin:
		if f[0] == '#': continue
		print f
		name, score, pacst = tuple( f.split()[:3] )
		name 	= name[0] + name[6] + name[7:12]
		score 	= float( score )
		pacst 	= float( pacst )

		if vacRMSD[name] < 1.8: 
			print name, '!!!!!!!!!!!!!!!!', vacRMSD[name]
			continue

		rmsd.append( vacRMSD[name] )
		a_sc.append( score )
		p_st.append( pacst )

rmsd = np.array( rmsd )
a_sc = np.array( a_sc )
p_st = np.array( p_st )

plt.plot( rmsd, a_sc, '.r' )
plt.ylabel('del_REU')
plt.xlabel('RMSD (Angstrom)') 
plt.title( 'CHAMP designs  vacRMSD (0.5 ns ) vs dimerREU\n'  )
plt.axis( [0, max( rmsd ) + 0.2, min( a_sc ) - 1, max( a_sc ) + 1] )
#path = os.path.join( dirPath, 'helix_rmsd.pdf' )
#plt.savefig( path, bbox_inches='tight' ) 
plt.show()

plt.plot( rmsd, p_st, '.r' )
plt.ylabel('del_REU')
plt.xlabel('RMSD (Angstrom)') 
plt.title( 'CHAMP designs  vacRMSD (0.5 ns ) vs dimerREU\n'  )
plt.axis( [0, max( rmsd ) + 0.2, min( p_st ) - 1, max( p_st ) + 1] )
#path = os.path.join( dirPath, 'helix_rmsd.pdf' )
#plt.savefig( path, bbox_inches='tight' ) 
plt.show()








