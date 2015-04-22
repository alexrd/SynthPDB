#!/usr/bin/env perl

use strict;

use lib "/Users/alexrd/perl";
use mj qw ( %names %symb );

my @excl = qw(mse yof);
our $cutdistsq = 36;

sub help {
    print("
listallligands.pl

Alex Dickson
University of Michigan

This script takes a folder full of pdbs, and lists all ligands, as well as a 
histogram of which residues they are closest to.  The PDBs should all have the 
same amino acid numbering (use alignbypattern.pl beforehand).

This program thus only looks for pdbs that have 'align' in the title.

It is assumed that all ligands will start with HETATM, and everything that starts with
HETATM is a ligand, except for entries in the 'excl' array, such as modified amino acids.

Current excl: (@excl)
");
    exit(0);
}

my $flag = shift(@ARGV);
if ($flag eq "-h" || $flag eq "-help") {
    &help;
}

my %ligands;  # $ligands{name}{closest} = count 
              # $ligands{name}{pdbs} = (pdb1, pdb2,..)

my @files = split(/\n/,`ls *align*.pdb`);

if (scalar @files == 0) {
    print("No aligned pdbs in this directory!\n");
    &help;
}

my %hetnam;
my %is_excl;
for (@excl) { $is_excl{$_} = 1; };

for my $file (@files) {

    my %prot;
    my %lig;

    open(MYIN,"$file");
    (my $alignchain) = ($file =~ m/align_([A-Z])/);
    (my $pdbname) = ($file =~ m/([a-zA-Z0-9]+)_align/);
    while (my $line = <MYIN>) {
	if (substr($line,0,4) eq "ATOM" || substr($line,0,6) eq "HETATM") {
	    my $chain = substr($line,21,1);
	    my $resnum = substr($line,22,4);
	    $resnum =~ s/^\s+//;
	    my $resname = substr($line,17,3);
	    $resname =~ s/^\s+//;	 

	    my @x = (substr($line,30,8),substr($line,38,8),substr($line,46,8));
	    for (@x) { $_ =~ s/^\s+//; }

	    if ($is_excl{lc($resname)} || substr($line,0,4) eq "ATOM") {
		$prot{$chain}{$resnum}{"name"} = lc($resname);
		$prot{$chain}{$resnum}{"pos"} = \@x;
#		print("adding res $resname with x = (@x) to prot chain $chain\n");
	    } else {
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
	    }	   
	} elsif (substr($line,0,6) eq "HETNAM") {
	    my $het = lc(substr($line,11,3));
	    $het =~ s/^\s+//;	 
	    my $name = substr($line,15,52);
	    $name =~ s/\s+$//;
	    if (!exists $hetnam{$het}{"names"}) {
		@{$hetnam{$het}{"names"}} = ();
	    }
	    if (&notin($hetnam{$het}{"names"},$name)) {
		push(@{$hetnam{$het}{"names"}},$name);
	    }
	}
	if (substr($line,0,6) eq "ENDMDL") {
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
		$ligands{$ligname}{"pdbs"}{$pdbname} = 1;
	    }
	}
    }
}

my $nlig = scalar( keys %ligands );

@{$hetnam{"hoh"}{"names"}} = qw(water);

print("A total of $nlig ligands were found:
-------------------------------------------\n");
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
	if ($res ne "tot" && $res ne "pdbs") {
	    print("  residue $res ($ligands{$l}{$res})\n");
	}
    }
    my @pdbs = sort keys %{$ligands{$l}{"pdbs"}};
    print("  PDBs:  @pdbs\n\n");
}       


sub getclosest {
    my @x = @{$_[0]};
    my %prot = %{$_[1]};

    my $mind = $cutdistsq;
    my $closeres = undef;
    for my $i ( keys %prot ) {
	my @x2 = @{$prot{$i}{"pos"}};
	my $dsq = 0;
	for my $j (0..@x-1) {
	    $dsq += ($x[$j] - $x2[$j])**2;
	}
	if ($dsq < $mind) {
	    $mind = $dsq;
	    $closeres = $i;
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
