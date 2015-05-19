#!/usr/bin/env perl

use strict;

use lib "/Users/alexrd/perl";
use mj qw ( %names %symb );

my @excl = qw(mse yof);
our $cutdistsq = 36;    # cutoff for residue-ligand COM interactions (squared)
our $cluster_rad = 5;
our $rescutoffsq = 16;  # cutoff for residue-residue interactions (squared)

my $alignfile = shift(@ARGV);
my $minres = shift(@ARGV);

sub help {
    print("
listallligands.pl

Usage:  ./listallligands.pl alignfile alignres

Alex Dickson
University of Michigan

This script takes a folder full of pdbs, and lists all ligands, as well as a 
histogram of which residues they are closest to.  The PDBs should all have the 
same amino acid numbering, and should be geometrically aligned.

The workflow is then:
1) groupalign.pl
2) alignstructures.pl
3) listallligands.pl

This program thus only looks for pdbs in the 'structalign/' directory.

The first argument is the alignment file used in groupalign.pl
The second argument is the index of alignment used in groupalign.pl.

It is assumed that all ligands will start with HETATM, and everything that starts with
HETATM is a ligand, except for entries in the 'excl' array, such as modified amino acids.
It uses the original pdbs to get the HETNAM entries, since these are changed to ATOM
by the VMD step in alignstructures.pl

Current excl: (@excl)

It uses a clustering radius to generate representative centers of geometry for each ligand,
which are saved as atoms in a pdb file, written to 'allligands.pdb'.

By default, this ignores water (which can swamp the results).  This can be included by 
uncommenting the line below that says 'uncomment to include water'.

For each aligned residue, it computes the fraction of contacts with residues outside of 
the domain.  These can either be on the same chain, or an a different protein altogether.
These are broken down into individual amino acids, and are printed in the table 'contact_freq.dat'
");
    exit(0);
}

my $flag = shift(@ARGV);
if ($flag eq "-h" || $flag eq "-help") {
    &help;
}

if (!-f $alignfile) {    
    print("Error:  alignfile does not exist!\n");
    &help;
}

my %alignseq;
# read alignfile
open(MYIN,$alignfile);
while (my $line = <MYIN>) {
    chomp($line);
    my @arr = split(/\s+/,$line);
    if (defined $arr[0] && $arr[0] ne "CLUSTAL") {
	$alignseq{$arr[0]} .= "$arr[1]";
    }
}
close(MYIN);

my @tkeys = keys %alignseq;
my $nres = length($alignseq{$tkeys[0]});
my $maxres = $minres + $nres - 1;  # residues from minres to maxres are in the aligned region

my @residues;
for my $i (0..$nres-1) {
    my %tmp;
    $residues[$i] = \%tmp;  # for each residue, keep a hash of what amino acids 
}                           # (outside the aligned region) it is closest to

my %ligands;  # $ligands{name}{closest_res} = count 
              # $ligands{name}{pdbs} = (pdb1, pdb2,..)
              # $ligands{name}{total} = total count
              # $ligands{name}{positions}[i][j] = x   (where i is the cluster index and j = 0,1,2 for x,y,z coord)

my @files = split(/\n/,`ls structalign/*.pdb`);

if (scalar @files == 0) {
    print("No pdbs in the structalign/ directory!\n");
    &help;
}

my %hetnam;
#@{$hetnam{"hoh"}{"names"}} = ("water");          # uncomment to include water
my %is_excl;
for (@excl) { $is_excl{$_} = 1; };

for my $file (@files) {

    my %prot;
    my %lig;

    # read through the original's HETNAM entries 
    (my $pdbname) = ($file =~ m/([a-zA-Z0-9]+)_struct_align/);
    my @lines = split(/\n/,`grep HETNAM ${pdbname}.pdb`);
    for my $line (@lines) {
	if (substr($line,0,6) eq "HETNAM") {
	    my $het = lc(substr($line,11,3));
	    $het =~ s/\s//g;	 
	    my $name = substr($line,15,52);
	    $name =~ s/\s//g;
	    if (!$is_excl{$het}) {
		if (!exists $hetnam{$het}{"names"}) {
		    @{$hetnam{$het}{"names"}} = ();
		}
		if (&notin($hetnam{$het}{"names"},$name)) {
		    push(@{$hetnam{$het}{"names"}},$name);
		}
	    }
	}
    }

    open(MYIN,"$file");
    (my $alignchain) = ($file =~ m/([A-Z])\.pdb/);
    while (my $line = <MYIN>) {
	if (substr($line,0,4) eq "ATOM" || substr($line,0,6) eq "HETATM") {
	    my $chain = substr($line,21,1);
	    my $resnum = substr($line,22,4);
	    $resnum =~ s/\s//g;
	    my $resname = substr($line,17,3);
	    $resname =~ s/\s//g;
	    my $atomname = substr($line,12,4);
	    $atomname =~ s/\s//g;
	    
	    my @x = (substr($line,30,8),substr($line,38,8),substr($line,46,8));
	    for (@x) { $_ =~ s/\s//g; }

	    if (exists $hetnam{lc($resname)}) {
		$lig{$chain}{$resnum}{"name"} = lc($resname);
		if (!exists $lig{$chain}{$resnum}{"pos"}) {
		    @{$lig{$chain}{$resnum}{"pos"}} = qw(0 0 0);
		    $lig{$chain}{$resnum}{"npos"} = 0;
		}
		for my $i (0..2) {
		    $lig{$chain}{$resnum}{"pos"}[$i] += $x[$i];
		}
		$lig{$chain}{$resnum}{"npos"}++;
#		print("adding res $resname with x = (@x) to lig chain $chain\n");
	    } else {
		$prot{$chain}{$resnum}{"name"} = lc($resname);
		$prot{$chain}{$resnum}{"pos"}{$atomname} = \@x;
#		print("adding res $resname with x = (@x) to prot chain $chain\n");
	    }
	} elsif (substr($line,0,6) eq "ENDMDL") {
	    last;
	}
    }
    close(MYIN);

    # for each ligand, find what prot residue is closest

    for my $c ( keys %lig ) {
	for my $n ( keys %{$lig{$c}} ) {
	    for (@{$lig{$c}{$n}{"pos"}}) {
		$_ /= $lig{$c}{$n}{"npos"};   # now this is the center of geometry of the ligand
	    }
	    my $xref = $lig{$c}{$n}{"pos"};
	    (my $n2) = &getclosest($xref,$prot{$alignchain});
	    if (defined $n2) {
		my $ligname = $lig{$c}{$n}{"name"};
		if (!exists $ligands{$ligname}{$n2}) {
		    $ligands{$ligname}{$n2} = 1;
		} else {
		    $ligands{$ligname}{$n2}++;
		}
		if (!exists $ligands{$ligname}{"tot"}) {
		    $ligands{$ligname}{"tot"} = 0;
		}
		$ligands{$ligname}{"tot"}++;
		if (!exists $ligands{$ligname}{"centers"}) {
		    my @tmp;
		    $ligands{$ligname}{"centers"} = \@tmp;
		}
		&getclusterind($xref,$ligands{$ligname}{"centers"},$pdbname);
	    }
	}
    }

    # for each residue in the aligned region, find out what it is touching, and save this to 
    #   the %residues hash

    for my $res ($minres..$maxres) {
	if (!exists $residues[$res-$minres]{"norm"}) {
	    $residues[$res-$minres]{"norm"} = 0;
	}		    
	$residues[$res-$minres]{"norm"}++;
	for my $c (keys %prot) {
	    for my $r (keys %{$prot{$c}} ) {
		if ($c ne $alignchain || ($r < $minres || $r > $maxres)) { # outside the aligned region
		    my $keepchecking = 1;
		    for my $at1 ( keys %{$prot{$alignchain}{$res}{"pos"}} ) {
			if ($keepchecking) {
			    for my $at2 ( keys %{$prot{$c}{$r}{"pos"}} ) {
				if ($keepchecking) {
				    my $dsq = &getdistsq($prot{$alignchain}{$res}{"pos"}{$at1},$prot{$c}{$r}{"pos"}{$at2});
				    if ($dsq < $rescutoffsq) {  # found a contact, record and move on
					if (!exists $residues[$res-$minres]{$prot{$c}{$r}{"name"}}) {
					    $residues[$res-$minres]{$prot{$c}{$r}{"name"}} = 0;
					}
					$residues[$res-$minres]{$prot{$c}{$r}{"name"}}++;
					$keepchecking = 0;
				    } elsif ($dsq > 400) {  # not even close, skip the rest of the atoms
					$keepchecking = 0;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

# print ligand results, write ligand coordinate centers to a pdb

my $nlig = scalar( keys %ligands );
print("A total of $nlig ligands were found:
-------------------------------------------\n");

open(MYOUT,">allligands.pdb");
my $atomcount = 1;

for my $l ( sort {$ligands{$b}{"tot"} <=> $ligands{$a}{"tot"}} keys %ligands) {
    print("ligand $l observed $ligands{$l}{'tot'} times\n");
    if (exists $hetnam{$l}{"names"}) {
	my @names = @{$hetnam{$l}{"names"}};
	print("  name(s): $names[0]\n");
	for my $i (1..@names-1) {
	    print("           $names[$i]\n");	
	}
    }
    for my $res ( sort {$ligands{$l}{$b} <=> $ligands{$l}{$a}} keys %{$ligands{$l}} ) {
	if ($res ne "tot" && $res ne "centers") {
	    print("  residue $res ($ligands{$l}{$res})\n");
	}
    }

    my $nclust = scalar @{$ligands{$l}{"centers"}};
    print("  Number of cluster centers: $nclust\n");
    for my $i (0..$nclust-1) {
	printf(MYOUT "ATOM    %3i CA   %3s A %3i    %8.3f%8.3f%8.3f  1.00  0.00           C\n",$atomcount,uc($l),$atomcount,$ligands{$l}{"centers"}[$i]{"pos"}[0],$ligands{$l}{"centers"}[$i]{"pos"}[1],$ligands{$l}{"centers"}[$i]{"pos"}[2]);
	printf("    Center $i: %f %f %f\n",$ligands{$l}{"centers"}[$i]{"pos"}[0],$ligands{$l}{"centers"}[$i]{"pos"}[1],$ligands{$l}{"centers"}[$i]{"pos"}[2]);
	printf("    PDBs: ");
	$atomcount++;
	my @pdbs = sort keys %{$ligands{$l}{"centers"}[$i]{"pdbs"}};
	for (@pdbs) {
	    printf("$_ ");
	}
	printf("\n");
    }
    print("\n");
}       
print(MYOUT "END\n");
close(MYOUT);

# print residue adjacency information to contact_freq.dat

open(MYOUT,">contact_freq.dat");
my %allkeys;
for my $i (0..@residues-1) {
    for (keys %{$residues[$i]}) {
	if ($_ ne "norm") {
	    $allkeys{$_} = 1;  # get all contacting residues
	}
    }
}
print(MYOUT "index ");  # print table column titles
for (sort keys %allkeys) {
    print(MYOUT "$_ ");
}
print(MYOUT "\n");
for my $res ($minres..$maxres) {
    my $i = $res - $minres;
    print(MYOUT "$res ");
    for (sort keys %allkeys) {
	if (!exists $residues[$i]{$_} || $residues[$i]{"norm"} == 0) {
	    printf(MYOUT "%f ", 0);
	} else {
	    printf(MYOUT "%f ", $residues[$i]{$_}/$residues[$i]{"norm"});
	}
    }
    print(MYOUT "\n");
}
close(MYOUT);

sub getclosest {
    my @x = @{$_[0]};
    my %prot = %{$_[1]};

    my $mind = $cutdistsq;
    my $closeres = undef;
    for my $i ( keys %prot ) {
	for my $at ( keys %{$prot{$i}{"pos"}} ) {
	    my @x2 = @{$prot{$i}{"pos"}{$at}};
	    my $dsq = 0;
	    for my $j (0..@x-1) {
		$dsq += ($x[$j] - $x2[$j])**2;
	    }
	    if ($dsq < $mind) {
		$mind = $dsq;
		$closeres = $i;
	    }
	}
    }
    
    return $closeres;
}
	    
sub notin {
    my @arr = @{$_[0]};
    my $el = $_[1];

    for (@arr) {
	if ($_ eq $el) {
	    return 0;
	}
    }
    return 1;
}

sub getclusterind {
    my @x = @{$_[0]};
    my $posarray = $_[1];
    my $pdb = $_[2];
    
    my $closest = undef;
    my $mind = $cluster_rad**2;
    my $n = @{$posarray};
    if ($n > 0) {
	for my $i (0..$n-1) {
	    my @x2 = @{$posarray->[$i]{"pos"}};
	    my $dsq = 0;
	    for my $j (0..@x-1) {
		$dsq += ($x[$j] - $x2[$j])**2;
	    }
	    if ($dsq < $mind) {
		$mind = $dsq;
		$closest = $i;
	    }
	}
    }

    if (!defined $closest) {  # define a new cluster center!
	my %tmp;
	my %pdbs = ($pdb => 1);
	$tmp{"pos"} = \@x;
	$tmp{"pdbs"} = \%pdbs;
	$posarray->[$n] = \%tmp;
    } else {
	$posarray->[$closest]{"pdbs"}{$pdb} = 1;
    }
}

sub getdistsq {
    my @x = @{$_[0]};
    my @y = @{$_[1]};

    my $d = 0;
    for my $i (0..@x-1) {
	$d += ($x[$i] - $y[$i])**2;
    }

    return $d;
}
