# pdb1 = with ligand
# pdb2 = no ligand
# source align_2structs.tcl
# pairrmsd 2o7n.pdb 1dgq.pdb "protein within 5 of ((not protein) and (not water)) and noh" rmsd.out

set i 0; foreach n $argv {set [incr i] $n}

set pdb1 $1
set pdb2 $2
set seltext1 $3
set seltext2 $4
set outfile $5
    
set seltext1 [string map {- " "} $seltext1]
set seltext2 [string map {- " "} $seltext2]

set output [open $outfile w]
    
mol new $pdb1 type pdb waitfor all
mol new $pdb2 type pdb waitfor all

set sel1 [atomselect 0 $seltext1 frame 0]
set sel2 [atomselect 1 $seltext2 frame 0]
set num1 [$sel1 num]
set num2 [$sel2 num]

if { $num1 == $num2 } then {
    set all [atomselect 0 all frame 0]

    set trans_mat [measure fit $sel1 $sel2]
    $all move $trans_mat
    set rmsd [measure rmsd $sel1 $sel2]
    puts $output "[format "%8.3f" $rmsd]"
    
    $all writepdb t.pdb 
} else {
    puts $output "[format "%8s" XXX]"

}
flush $output
close $output

