#!/usr/bin/python

#$ -S /usr/bin/python 
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=00:05:00
#$ -cwd
#$ -j y
#$ -o /netapp/home/mmravic/CHAMP/EGFR/bbDesign_CTERMrnd2/logfiles 
#$ -t 1-10000


# Marco Mravic DeGrado Lab 2016


# input 1: pdb to input, oriented in the membrane.. auto grabs resfile and span file  
# input 2: protocol xml file
# input 3: path to rosetta 

# output: directory called 'designs' to store output files for design trajactories

## example command line
# python ~/CHAMP/bin/champDesign_rosMem.py ~/CHAMP/a5/bb2Design/match_1/match_1.pdb ~/CHAMP/bin/helix_Design.xml ~/rosetta/



import sys, os, subprocess as sp, time
start = time.time()

# I/O

inPdb 		= sys.argv[1]
inSpan		= inPdb[:-4] + '.span'
inXML		= sys.argv[2]
inResF		= os.path.join( os.path.dirname(inPdb), 'resfile' )

rosiBase 	= sys.argv[3]
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiDB 		= os.path.join( rosiBase, 'database/' )

oDir 		= os.path.join( os.path.dirname(inPdb), 'outputs/' )

try:
	output_suffix 			= '_out%s' % (str(  os.environ["SGE_TASK_ID"]) )
except KeyError:
	output_suffix                      = '_out%s' % ( 'Local' )

if not os.path.exists( oDir ):
	os.mkdir( oDir )

#

n = '1'

cmd = [  rosiScrps, 
'-parser:protocol', inXML, 					# Path to Rosetta script (see above)
'-in:file:s', inPdb,						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:setup:spanfiles', inSpan,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix', oDir,
'-out:suffix', output_suffix,
'-packing:resfile', inResF,
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

print 
print cmd
print 

sp.call( cmd )

print 
print 'Entire run took this many seconds:', time.time() - start
