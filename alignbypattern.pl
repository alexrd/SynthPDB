#!/usr/bin/env perl

use strict;

use lib "/Users/alexrd/perl";
use mj qw ( %names %symb );

my @residue = split(/:/,shift(@ARGV));
my @index = split(/:/,shift(@ARGV));

sub help {
    print("
alignbypattern.pl

Alex Dickson
University of Michigan

This script takes a folder full of pdbs, and aligns them according to a pattern of conserved residues.
If the pattern is found more than once it will align to each instance, and then print a warning.

The first argument is a colon-separated list of the residues that are conserved using one-letter codes.
Multiple residues can be specified to allow for flexibility:
e.g. Y:N:ILVM:DN 

The second argument is a colon-separated list of the corresponding residue indices.
e.g. 412:428:425:421
");
    exit(0);
}

if ((scalar @residue != scalar @index) || (scalar @residue == 0)) {    	  
    print("Usage:  res1:res2:res3:.. ind1:ind2:ind3:..\n");
    &help;
}

my @files = split(/\n/,`ls *.pdb`);

if (scalar @files == 0) {
    print("No pdbs in this directory!\n");
    &help;
}

for my $file (@files) {

    my %resi;

    open(MYIN,"$file");
    while (my $line = <MYIN>) {
	if (substr($line,0,4) eq "ATOM" || substr($line,0,6) eq "HETATM") {
	    my $chain = substr($line,21,1);
	    my $resnum = substr($line,22,4);
	    $resnum =~ s/^\s+//;
	    my $resname = substr($line,17,3);
	    $resi{$chain}{$resnum} = lc($resname);
	}
	if (substr($line,0,6) eq "ENDMDL") {
	    last;
	}
    }
    close(MYIN);

    my @matchchain = ();
    my @matchnum = ();
    for my $c ( keys %resi ) {
	for my $n ( keys %{$resi{$c}} ) {
	    if (exists $symb{$resi{$c}{$n}}) {
		my $aa = $symb{$resi{$c}{$n}};
		if ($residue[0] =~ m/$aa/) {   # the first one matches, check the rest

#		    print("following up on chain $c num $n : resid = $aa\n");

		    my $match = 1;
		    for my $i (1..@residue-1) {
			if ($match) {
			    my $offset = $index[$i] - $index[0];
			    my $ni = $n + $offset;
			    if (exists $resi{$c}{$ni}) {
				my $aai = $symb{$resi{$c}{$ni}};
				if (!($residue[$i] =~ m/$aai/)) {
#				    print("does not match because of resid $ni : $resi{$c}{$ni} : ($aai instead of $residue[$i])\n");
				    $match = 0;
				}
			    } else {
#				print("resi $resi{$c}{$ni} does not exist in the table\n");
				$match = 0;
			    }
			}		   
		    }
		    if ($match) {
			push(@matchchain,$c);
			push(@matchnum,$n);
		    }
		}
	    }
	}
    }

    if (scalar @matchchain == 0) {
	print(STDERR "no sequence matches found for $file\n");
    } else {
	my $senterr = 0;
	my %chains;
	for my $i (0..@matchchain-1) {
	    if (!exists $chains{$matchchain[$i]}) {
		$chains{$matchchain[$i]} = 1;
	    } else {
		$chains{$matchchain[$i]}++;
		if (!$senterr) {
		    print(STDERR "multiple sequence matches found for $file ");
		    for my $j (0..@matchchain-1) {
			print(STDERR "($matchchain[$j])$matchnum[$j] ");
		    }
		    print(STDERR "\n");
		    $senterr = 1;
		}
	    }
	    my $outfile = $file;
	    $outfile =~ s/\.pdb/_align_${matchchain[$i]}$chains{$matchchain[$i]}\.pdb/;

	    open(MYOUT,">$outfile");
	    open(MYIN,"$file");
	    while (my $line = <MYIN>) {
		if ((substr($line,0,4) eq "ATOM" || substr($line,0,4) eq "TER ") || substr($line,0,6) eq "HETATM") { 
		    my $chain = substr($line,21,1);
		    my $resnum = substr($line,22,4);
		    if ($chain eq $matchchain[$i]) {
			my $newnum = $resnum - $matchnum[$i] + $index[0];
			if ($newnum < 0) {
			    $newnum = $resnum;
			    substr($line,21,1,"X");
			} else {
			    my $strnum = sprintf('%4s',$newnum);
			    substr($line,22,4,$strnum);
			}
			print(MYOUT $line);
		    } else {
			print(MYOUT $line);
		    }
		} else {
		    print(MYOUT $line);
		}
	    }
	    close(MYIN);
	    close(MYOUT);
	}
    }
}
