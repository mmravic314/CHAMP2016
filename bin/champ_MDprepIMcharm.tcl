#champ_MDprep.tcl
package require autopsf
package require ssrestraints

set id [mol new {psf_Input.pdb}]
autopsf -mol $id -prefix champ -protein -top /Users/mmravic/bin/toppar/top_all36_prot.rtf

ssrestraints -psf champ.psf -pdb champ.pdb -o dihedrals.cst -hbonds -k_prot 300 -hbbondk 30
exit