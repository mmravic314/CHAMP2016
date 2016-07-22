# python ~/CHAMP/bin/prep4MD.py ~/CHAMP/a5/model_list.txt ~/CHAMP/a5/MD_runs/ ~/CHAMP/bin/vacMD.conf ~/CHAMP/bin/champ_MDprep.tcl ~/CHAMP/bin/Z_align.pdb ~/bin/NAMD/namd2 /Applications/VMD1.9.2.app/Contents/MacOS/startup.command ~/CHAMP/bin/vacEQ.conf 

# transfer a list of CHAMP files
 # prep as psf & MD input
import sys, os, shutil, subprocess as sp, numpy as np, time
from prody import *

# I/O
pdbList = open( sys.argv[1], 'rU' ).readlines()
targDir = sys.argv[2]			# path to directory that MD files will be generated and output into

conFmd = sys.argv[3]			# path to generic MD run configuration file, no constraints
tcl_prep= sys.argv[4]			# path to tcl script for VMD to generate psf/pdb inputs for MD runs

alignZ	= parsePDB( sys.argv[5] ) 	# path to pesudoatoms oriented at "bilayer" z = -15,0,15 & {x,y} = -5,0,5
namdPath= sys.argv[6]			# path to namd2 binary
vmdPath	= sys.argv[7]			# path to VMD startup.command (binary)
conFEQ	= sys.argv[8]			# path to generic minimization & equilibration config file, constraints

FNULL= open(os.devnull, 'w')	# path to dump verbose standard output to
#

for f in pdbList:
	if f[0] == '#': continue

	## Transfer files and make working directory
	spl = os.path.basename(f.rstrip())[:-4].split('_')
	modelID 	= '%s%s_%s' % ( spl[0][0], spl[1], spl[2] ) 
	targInner 	= os.path.join( targDir, modelID )
	
	#print targInner, spl, f
	#sys.exit()

	if not os.path.exists( targInner ): 
		os.mkdir( targInner )

	# Skip if simulation alread looks complete?

	## code here for that, if you want...


	# remove membrane residue, rename, and transform protein to membrane normal using MEM residue (ROSETTA object)
	inPDB 				= parsePDB( f.rstrip() )
	# extract center and normal vetor to add pseudo atoms +/-15 from center...
	mem 				= inPDB.select( 'resname MEM' )

	# If not a rosetta output file, assume it is already oriented into a bilayer slab
	if not mem:
		# Gen PSF file for dimer


		os.chdir( targInner )
		writePDB( os.path.join( targInner, 'psf_Input.pdb'), inPDB )

		cmd = [ vmdPath, '-dispdev', 'text', '-e', tcl_prep ]
		sp.Popen( cmd , stdout=FNULL)
		print 'Input files created for', targInner

		# prepare constraint file
		tmp = parsePDB('champ.pdb')
		for k in tmp.select('ca'):	# original file should be usable as constraint file with CA restained initially
			k.setBeta( 1.0 )
		writePDB( 'champ.pdb', tmp)

		shutil.copy( conFEQ, targInner )
		shutil.copy( conFmd, targInner )

		md_start = time.time()
		## Run simulations
		eqCMD = [ namdPath, os.path.basename( conFEQ ) ]		# print out to debug
		prCMD = [ namdPath, os.path.basename( conFmd ) ]		# print out to debug

		eqLog = sp.Popen( eqCMD , stdout=sp.PIPE).communicate()[0]
		eqOut = open( 'outputEQ.log', 'w' )
		eqOut.write( eqLog)
		del eqLog
		eqOut.close()
		print 'EQ run complete for', targInner

		prLog = sp.Popen( prCMD , stdout=sp.PIPE).communicate()[0]
		prOut = open( 'outputPROD.log', 'w' )
		prOut.write( prLog )
		del prLog
		prOut.close()
		print 'Production run complete for', targInner, '\ntime taken (s):', time.time()-md_start, '\n'

		continue

	thick, cent, norm 	= tuple( mem.getCoords() )
	top, bottom			= cent + 15*( norm - cent ), cent - 15*( norm - cent )
	topX, bottomX		= cent + 15*( norm - cent ) + [5,0,0], cent - 15*( norm - cent ) + [5,0,0]
	topY, bottomY		= cent + 15*( norm - cent ) + [0,5,0], cent - 15*( norm - cent ) + [0,5,0]

	norm2 = norm - cent

	marks 				= alignZ.copy()
	cSet = np.array( [top, cent, bottom, topX, bottomX, topY, bottomY] )
	marks.setCoords( cSet )


	# Add pseudo atoms to object, calculate transformation matrix to align these to [(0,0,15), (0,0,0), (0,0,-15)]
	# Then write to file without these pseudo atoms
	inPDB += marks
	marks = inPDB.select( 'name ZN' )
	trans = calcTransformation( marks.getCoords(), alignZ )
	applyTransformation( trans, inPDB )
	inPDB = inPDB.select( '(not name ZN) (not resname MEM) heavy' ).copy()
	inPDB.setTitle( f.rstrip() )
	writePDB( os.path.join( targInner, 'psf_Input.pdb'), inPDB )


	# Gen PSF file for dimer
	os.chdir( targInner )

	cmd = [ vmdPath, '-dispdev', 'text', '-e', tcl_prep ]
	sp.Popen( cmd , stdout=FNULL)
	print 'Input files created for', targInner

	# prepare constraint file
	tmp = parsePDB('champ.pdb')
	for k in tmp.select('ca'):	# original file should be usable as constraint file with CA restained initially
		k.setBeta( 1.0 )
	writePDB( 'champ.pdb', tmp)

	shutil.copy( conFEQ, targInner )
	shutil.copy( conFmd, targInner )


	## Run simulations
	md_start = time.time()

	eqCMD = [ namdPath, os.path.basename( conFEQ ) ]		# print out to debug
	prCMD = [ namdPath, os.path.basename( conFmd ) ]		# print out to debug

	eqLog = sp.Popen( eqCMD , stdout=sp.PIPE).communicate()[0]
	eqOut = open( 'outputEQ.log', 'w' )
	eqOut.write( eqLog)
	del eqLog
	eqOut.close()
	print 'EQ run complete for', targInner

	prLog = sp.Popen( prCMD , stdout=sp.PIPE).communicate()[0]
	prOut = open( 'outputPROD.log', 'w' )
	prOut.write( prLog )
	del prLog
	prOut.close()
	print 'Production run complete for', targInner, '\ntime taken (s):', time.time()-md_start, '\n'

	#sys.exit()


''' 
sp.call( [ 'vmd', '-dispdev', 'text', '-e', setUp ] )

>Main< (mmravic) 21 % autopsf
Usage: autopsf -mol <molnumber> <option1> <option2>...
Options:
  Switches:
    -protein (only generate psf for protein)
    -nucleic (only generate psf for nucleic acid segment)
    -solvate (run solvate on structure after psf generation)
    -ionize (run autoionize on structure after psf generation)
    -noguess (don't use guesscoord during psf building)
    -noterm (don't automatically add terminii to proteins)
    -rotinc <increment> (degree increment for rotation)
    -gui (force graphical interface mode)
    -regen (regenerate angles/dihedrals)
    -splitonly (split chains, but don't build structure)
  Single option arguments:
    -mol <molecule> (REQUIRED: Run on molid <molecule>)
    -prefix <prefix> (Use <prefix> as the prefix to output files)
    -top <top> (Use topology file <top>)
    -include <selection> (Include fragment specified by <selection> in psf generation
    -patch <patch> (Apply a patch -- see docs for syntax)

HETATM  978 THKN MEM Y  59      19.484   2.532   0.539  1.00  0.00           X  
HETATM  979 CNTR MEM Y  59       6.439  -4.852  -0.032  1.00  0.00           X  
HETATM  980 NORM MEM Y  59       6.415  -4.885   0.967  1.00  0.00           X 

'''