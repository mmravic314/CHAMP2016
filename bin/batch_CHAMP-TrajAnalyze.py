import sys, os, numpy as np, subprocess as sp, time
from prody import *

import matplotlib.pyplot as plt


# In Vacuum runs
# python ~/CHAMP/bin/batch_CHAMP-TrajAnalyze.py ~/CHAMP/a5/MD_runs/ /Applications/VMD1.9.2.app/Contents/MacOS/startup.command ~/CHAMP/bin/sasa.tcl

# IMM1 runs, go to ca. line 130 to change some options for analysis
# python ~/CHAMP/bin/batch_CHAMP-TrajAnalyze.py /Volumes/DJ\ MCHERRY/CHARMM_GBSW/gbsw_dynamics/ /Applications/VMD1.9.2.app/Contents/MacOS/startup.command ~/CHAMP/bin/sasa_IMM1.tcl


vmdPath	= sys.argv[2]			# path to VMD startup.command (binary)
tcl_sasa= sys.argv[3]
FNULL= open(os.devnull, 'w')	


dataStore = {}
for sdir in os.listdir( sys.argv[1] ):


	

	dirPath = os.path.join( sys.argv[1], sdir  )

	if not os.path.isdir( dirPath ): continue
	os.chdir( dirPath )

	def IMM1_parse( dirPath ):

		dcd = os.path.join( dirPath, 'dimer/gbswmemb1_01.dcd' )
		pdb = os.path.join( dirPath, 'dimer/gbswmemb1_01.pdb' ) 	# final frame PDB


		print pdb
		print dcd

		ensemble 	= parseDCD( dcd )
		pdb 		= parsePDB( pdb )


		# for imm1, reset the pdb object to the fist frame
		ensemble.setAtoms( pdb )
	#	ensemble.setCoords( pdb ) 
		pdb		= ensemble[1].getAtoms()



		## 
		# generically grab everything but 4 terminal residues of each helix
		h1, h2 	= set( pdb.getSegnames() )
		h1_resi = ' '.join( [ str(x) for x in  pdb.select( 'ca segment %s' % h1 ).getResnums()[4:-4] ] ) 
		h2_resi = ' '.join( [ str(x) for x in  pdb.select( 'ca segment %s' % h2 ).getResnums()[4:-4] ] ) 
		selStr 	= 'ca resnum %s %s ' % ( h1_resi, h2_resi )
		intface = pdb.select( selStr )
	
		# calc RMSD of interface
		ensemble.setAtoms(  intface )
		#ensemble.setCoords( intface ) 


		ensemble.superpose()
		rmsd = ensemble.getRMSDs()[:]

		## Write RMSD data to txt file
		stp 	= 0
		txt 	= ''
		for k in rmsd:
			txt += '%s\t%.3f\n' % ( stp, k )
			stp +=1
		oF 		= open( 'rmsd.dat', 'w' )
		oF.write( txt )
		oF.close() 


		## 		

		# call vmd to get SASA, path to script (if it isnt already done yet)

		record = time.time()

		print 'working on sata.dat for %s' % dirPath,

		if not os.path.exists('sasa.dat') or os.path.getsize( 'sasa.dat' ) < 50:
			cmd = [ vmdPath, '-dispdev', 'text', '-e', tcl_sasa ]
			sp.Popen( cmd , stdout=FNULL)
			time.sleep(3)

			while os.path.getsize( 'sasa.dat' ) < 50 or not os.path.exists('sasa.dat'):
				time.sleep(1)

		print os.path.getsize( 'sasa.dat' )

		record = time.time() - record
		print '... took %0.3f s' % record


		# load & analyze cat
		dSASA = []
		with open( 'sasa.dat' ) as fin:
			for f in fin:
				dSASA.append( round( float( f.rstrip() ), 3 ) )

		dSASA = np.array( dSASA )
		#print len(  dSASA )
		# store the mean and stdev change in SASA of dimerization: first 5 frames and final 50 frames
		dSASA_init_mean 	= np.mean( 	dSASA[:15] )
		dSASA_init_stdev 	= np.std( dSASA[:15] )
		dSASA_final_mean 	= np.mean( 	dSASA[-50:] )
		dSASA_final_stdev	= np.std( dSASA[-50:] )



		# plot RMSD
		plt.plot( np.arange( 1, len(rmsd) + 1 ), rmsd, '.r-' )
		plt.xlabel('Frame Number, (2.5 ps)')
		plt.ylabel('RMSD (Angstrom)') 
		plt.title( 'CHAMP designs C-alpha RMSD, 1 ns\n%s' %  sdir  )
		plt.axis( [0, 200, 0, 4.5] )
		path = os.path.join( dirPath, 'helix_rmsd.pdf' )
		plt.savefig( path, bbox_inches='tight' ) 
		#plt.show()

		plt.close()

		# decide if the simulation converged
		mean, stdev = np.mean( rmsd[20:] ), np.std( rmsd[20:] )
		# filePath RMSD_mean RMSD_stdev dSASA_init_mean dSASA_init_stdev dSASA_final_mean dSASA_final_stdev
		data = '%s %0.3f  %0.3f %0.3f  %0.3f %0.3f  %0.3f' % (sdir, mean, stdev, dSASA_init_mean, dSASA_init_stdev, dSASA_final_mean, dSASA_final_stdev) 

		dataStore[ sdir ] = data


		os.chdir( sys.argv[1] )
		print data, '\n\n'


		return


	### HASH this option on for IMM1 simulations
	IMM1_parse( dirPath )

	continue

	
#	dcd = os.path.join( dirPath, 'output_prod.dcd' )
#	pdb	= os.path.join( dirPath, 'champ.pdb' )
#	log	= os.path.join( dirPath, 'outputPROD.log' )
	
	ensemble 	= parseDCD( dcd )
	pdb 		= parsePDB( pdb )


	# for imm1, reset the pdb object to the fist frame
	ensemble.setAtoms( pdb )
	pdb		= ensemble[0].getAtoms()
	ensemble.setCoords( pdb ) 


	## 
	# generically grab everything but 4 terminal residues of each helix
	h1, h2 	= set( pdb.getChids() )
	h1_resi = ' '.join( [ str(x) for x in  pdb.select( 'ca chain %s' % h1 ).getResnums()[4:-4] ] ) 
	h2_resi = ' '.join( [ str(x) for x in  pdb.select( 'ca chain %s' % h2 ).getResnums()[4:-4] ] ) 
	selStr 	= 'ca resnum %s %s ' % ( h1_resi, h2_resi )
	intface = pdb.select( selStr )
	
	# calc RMSD of interface
	ensemble.setAtoms( intface )
	ensemble.superpose()
	rmsd = ensemble.getRMSDs()[-100:]

	# call vmd to get SASA, path to script
	cmd = [ vmdPath, '-dispdev', 'text', '-e', tcl_sasa ]
	#print ' '.join( cmd )
	#sp.Popen( cmd , stdout=FNULL)

	#time.sleep(50)

	# load & analyze cat
	dSASA = []
	with open( 'sasa.dat' ) as fin:
		for f in fin:
			dSASA.append( round( float( f.rstrip() ), 3 ) )

	dSASA = np.array( dSASA )
	#print len(  dSASA )
	# store the mean and stdev change in SASA of dimerization: first 5 frames and final 50 frames
	dSASA_init_mean 	= np.mean( 	dSASA[:15] )
	dSASA_init_stdev 	= np.std( dSASA[:15] )
	dSASA_final_mean 	= np.mean( 	dSASA[-50:] )
	dSASA_final_stdev	= np.std( dSASA[-50:] )



	# plot RMSD
	plt.plot( np.arange( 1, len(rmsd) + 1 ), rmsd, '.r-' )
	plt.xlabel('Frame Number, (2.5 ps)')
	plt.ylabel('RMSD (Angstrom)') 
	plt.title( 'CHAMP designs C-alpha RMSD, 1 ns\n%s' %  sdir[:-1]  )
	plt.axis( [0, 200, 0, 4.5] )
	path = os.path.join( dirPath, 'helix_rmsd.pdf' )
	plt.savefig( path, bbox_inches='tight' ) 
#	plt.show()

	plt.close()

	# decide if the simulation converged
	mean, stdev = np.mean( rmsd[20:] ), np.std( rmsd[20:] )
	# filePath RMSD_mean RMSD_stdev dSASA_init_mean dSASA_init_stdev dSASA_final_mean dSASA_final_stdev
	data = '%s %0.3f  %0.3f %0.3f  %0.3f %0.3f  %0.3f' % (sdir, mean, stdev, dSASA_init_mean, dSASA_init_stdev, dSASA_final_mean, dSASA_final_stdev) 
	#print data

	dataStore[ sdir ] = data


	os.chdir( sys.argv[1] )
	print data, '\n\n'

print 'filePath RMSD_mean RMSD_stdev dSASA_init_mean dSASA_init_stdev dSASA_final_mean dSASA_final_stdev [n=200]'
for k, v in dataStore.items():
	
	print v

