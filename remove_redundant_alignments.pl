#!/usr/bin/env perl

my $alignfile = shift(@ARGV);

sub help {
    print("
remove_redundant_alignments.pl

Usage:  remove_redundant_alignments.pl alignfile > newalignfile

Alex Dickson
University of Michigan

This script reads a set of alignments and removes the redundant ones.  It also removes sequences that have 
more general cases already defined.

The first argument is the name of the file where the alignment is stored.  The format of this file is assumed to
be the same as Clustal Omega output, which is as follows:

Q6IRX1_68       -----------------------LPDYY------LTIKKPMDMEKIRSHMMA---NKYQDIDSMVEDFVMMFNNACTYNEP--------
P25440_317      ----KHAAYAWPFYKPVDASALGLHDYH------DIIKHPMDLSTVKRKMEN---RDYRDAQEFAADVRLMFSNCYKYNPPDHDVV---
Q58F21_287      ----KHFSYAWPFYNPVDVNALGLHNYY------DVVKNPMDLGTIKEKMDN---QEYKDAYKFAADVRLMFMNCYKYNPPDHEVV---
B4E3L4_723      -----------PFRQPVDPQLLGIPDYF------DIVKNPMDLSTIKRKLDT---GQYQEPWQYVDDVWLMFNNAWLYNR---------
P51532_1461     -----------------------LPEYY------ELIRKPVDFKKIKERIRN---HKY-------------------------------
etc..

Where the first word is a unique identifier (here, the uniprot ID followed by the starting position)
and the second word is the alignment for each sequence.  Multiple 'paragraphs' can be used for longer sequences, and 
they should be separated by two blank lines.  The first word identifier in each line of the second paragraph should 
be the same as was given in the first paragraph:  the start residue does not need to be updated.

");
    exit(0);
}

if (!-f $alignfile) {    
    print("Error:  alignfile does not exist!\n");
    &help;
}

my %alignseq;

open(MYIN, "$alignfile");
while (my $line = <MYIN>) {
    chomp($line);
    my @arr = split(/\s+/,$line);
    if (defined $arr[0] && $arr[0] ne "CLUSTAL") {
	$alignseq{$arr[0]} .= "$arr[1]";
    }
}
close(MYIN);

# remove exact entries
my %a2;
for my $u1 ( keys %alignseq ) {
    if (exists $a2{$alignseq{$u1}}) {
	delete $alignseq{$u1};
    } else {
	$a2{$alignseq{$u1}} = $u1;
    }
}

my %nosubsets = %alignseq;

# remove subsets
for my $u1 (keys %alignseq) {
    my @s1 = split(//,$alignseq{$u1});
    for my $u2 (keys %alignseq) {
	my $nd1 = ($alignseq{$u1} =~ m/(-)/g);
	my $nd2 = ($alignseq{$u2} =~ m/(-)/g);
	if ($nd1 > $nd2) {  # u1 has more dashes:  is it a subset of u2?
	    my @s2 = split(//,$alignseq{$u2});
	    my $match = 1;
	    for my $i (0..@s1-1) {
		if ($s1[$i] ne "-" && $s1[$i] ne $s2[$i]) {
		    $match = 0;
		    last;
		}
	    }
	    if ($match) {  # u1 is a subset, delete u2
		delete $nosubsets{$u2};
	    }
	}
    }
}

for ( keys %nosubsets ) {
    print("$_\t$nosubsets{$_}\n");
}
    
