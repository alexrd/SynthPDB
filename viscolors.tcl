# pdb1 = with ligand
# pdb2 = no ligand
# source align_2structs.tcl
# pairrmsd 2o7n.pdb 1dgq.pdb "protein within 5 of ((not protein) and (not water)) and noh" rmsd.out

set i 0; foreach n $argv {set [incr i] $n}

set phobpdb $1
set philpdb $2
set protpdb $3
set waterpdb $4
set exactpdb $5
set typepdb $6
set firstres $7
set lastres $8
set chain $9

mol new "allligands.pdb" type pdb waitfor all
mol new $phobpdb type pdb waitfor all
mol new $philpdb type pdb waitfor all
mol new $protpdb type pdb waitfor all
mol new $waterpdb type pdb waitfor all
mol new $exactpdb type pdb waitfor all
mol new $typepdb type pdb waitfor all

mol modstyle 0 0 VDW 0.400000 12.000000
mol modcolor 0 0 ResName

mol rename 1 {hydrophobic contacts}
mol rename 2 {hydrophilic contacts}
mol rename 3 {total protein contacts}
mol rename 4 {water contacts}
mol rename 5 {seq. cons. exact}
mol rename 6 {seq. cons. type}

for {set i 1} {$i < 7} {incr i} {
    mol modselect 0 $i resid $firstres to $lastres and chain $chain
    mol modcolor 0 $i Occupancy
    mol modstyle 0 $i NewCartoon 0.300000 10.000000 4.100000 0
    mol addrep $i
    mol modstyle 1 $i CPK 1.000000 0.300000 10.000000 10.000000
    mol modcolor 1 $i Occupancy
    mol modselect 1 $i resid $firstres to $lastres and chain $chain
    mol addrep $i
    mol modcolor 2 $i ColorID 2
    mol modstyle 2 $i NewCartoon 0.300000 10.000000 4.100000 0
    mol modselect 2 $i protein and not (resid $firstres to $lastres and chain $chain)
    if { $i == 1 } then {
	mol on $i
    } else {
	mol off $i
    }
}

