#sasa.tcl
# /Applications/VMD1.9.2.app/Contents/MacOS/startup.command -dispdev text -e sasa.tcl
set id [mol new output_prod.dcd waitfor all]
mol addfile champ.pdb

set group1 [atomselect $id "resid 7 to 26"]
set group2 [atomselect $id "resid 36 to 54"]
set group3 [atomselect $id "resid 7 to 26 36 to 54"]

set n_frames [molinfo top get numframes]

set f [open sasa.dat w]

for {set i 0} {$i < $n_frames} {incr i} {
	$group1 frame $i
	$group2 frame $i
	$group3 frame $i
	set sasa1 [measure sasa 1.4 $group1]
	set sasa2 [measure sasa 1.4 $group2]
	set sasa3 [measure sasa 1.4 $group3]
	set sasa [expr  $sasa2  + $sasa1 - $sasa3]
	puts $f $sasa
}

close $f

$group1 delete
$group2 delete
$group3 delete

exit

# measure sasa 1.4 $group