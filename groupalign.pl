#!/usr/bin/env perl

use strict;

use lib "/Users/alexrd/perl";
use mj qw ( %names %symb );

my $alignfile = shift(@ARGV);
my $alignres = shift(@ARGV);

sub help {
    print("
groupalign.pl

Usage:  groupalign.pl alignfile alignres

Alex Dickson
University of Michigan

This script uses a multiple-sequence alignment, and maps the indices of amino acids in PDBs such that they agree with 
this alignment.  It can handle alignments of domains (such as bromodomains, zinc fingers, etc.) where multiple 
domains can occur on the same protein.  The output pdbs are labeled by which chain was matched to the given sequence.

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

The second argument is the index if the first amino acid in the alignment.  This will be shared among all the 
aligned PDB files.  If during the alignment process, negative indices are obtained, these residues are put on 
a new chain, X, and their residue numbers are left unchanged.  To avoid this, consider using a higher index value.
");
    exit(0);
}

if (!-f $alignfile) {    
    print("Error:  alignfile does not exist!\n");
    &help;
}

my @files = split(/\n/,`ls *.pdb`);

if (scalar @files == 0) {
    print("No pdbs in this directory!\n");
    &help;
}

if (!-d "seqalign") {
    mkdir("seqalign");
}

# read alignment file

my %alignseq;
my %ndashbef;  # stores the number of dashes before each amino acid

open(MYIN,$alignfile);
while (my $line = <MYIN>) {
    chomp($line);
    my @arr = split(/\s+/,$line);
    if (defined $arr[0] && $arr[0] ne "CLUSTAL") {
	$alignseq{$arr[0]} .= "$arr[1]";
    }
}
close(MYIN);

my %seq = %alignseq;
for my $id ( keys %seq ) {
    $seq{$id} =~ s/-//g;  # get sequences with no dashes
    
    my @tmp;
    my $n = length($alignseq{$id});
    my $nres = 0;
    my $nd = 0;  # running number of dashes
    my $naa = 0;  # running number of amino acids
    
    for my $i (0..$n-1) {
	if (substr($alignseq{$id},$i,1) eq "-") {
	    $nd++;
	} else {
	    $tmp[$naa] = $nd;
	    $naa++;
	}
    }
    $tmp[$naa] = $nd;
    $ndashbef{$id} = \@tmp;
}

my $shiftupnexttime = 0;
my %flagsent;

for my $file (@files) {

    my %resi;

    my %lowres; # lowest residue index on each chain
    my %highres; # highest residue index on each chain

    open(MYIN,"$file");
    while (my $line = <MYIN>) {
	if (substr($line,0,4) eq "ATOM" || substr($line,0,6) eq "HETATM") {
	    my $chain = substr($line,21,1);
	    my $resnum = substr($line,22,4);
	    $resnum =~ s/\s+//;

	    if (!defined $lowres{$chain} || $resnum < $lowres{$chain}) {
		$lowres{$chain} = $resnum;
	    }
	    if (!defined $highres{$chain} || $resnum > $highres{$chain}) {
		$highres{$chain} = $resnum;
	    }

	    my $resname = substr($line,17,3);
	    $resi{$chain}{$resnum} = lc($resname);
	}
	if (substr($line,0,6) eq "ENDMDL") {
	    last;
	}
    }
    close(MYIN);

    # try to match each chain with sequence fragments
    my @match = ();
    for my $c ( keys %resi ) {

	# make a sequence "word" with spaces for missing residues
	$resi{$c}{"seq"} = "";
	for my $res ($lowres{$c}..$highres{$c}) {
	    my $resname = $resi{$c}{$res};
	    if (defined $resname) {
		my $aa = $symb{lc($resname)};
		if (!defined $aa) {
		    $aa = "X";
		    if (!exists $flagsent{$resname}) {
#			print("No residue entry for residue $resname\n");
		    }
		    $flagsent{$resname} = 1;
		}
		$resi{$c}{"seq"} .= $aa;		
	    } else {
		$resi{$c}{"seq"} .= " ";
	    }
	}

	# look through the list of sequence fragments for a match
	for my $s ( keys %seq ) {    
	    my $mref = &getmatch($resi{$c}{"seq"},$seq{$s});
	    my @tmatch = @{$mref};
	    if (scalar @tmatch > 0) {
		my $nm = scalar @tmatch;
		foreach (@tmatch) {
		    $_->{"ID"} = $s;
		    $_->{"chain"} = $c;
		    $_->{"startres"} += $lowres{$c};
		    push(@match,$_);  # collect matches in @match array
		}
	    }
	    
	}
    }

    if (scalar @match == 0) {
	print(STDERR "no sequence matches found for $file\n");
    } else {
	my $senterr = 0;
	my %chains;
	for my $i (0..@match-1) {
	    if (!exists $chains{$match[$i]{"chain"}}) {
		$chains{$match[$i]{"chain"}} = 1;
	    } else {
		$chains{$match[$i]{"chain"}}++;
		if (!$senterr) {
		    print(STDERR "multiple sequence matches found for $file ");
		    for my $j (0..@match-1) {
			print(STDERR "$match[$j]{'ID'}:($match[$j]{'chain'})$match[$j]{'startres'}-$match[$j]{'strength'} ");
		    }
		    print(STDERR "\n");
		    $senterr = 1;
		}
	    }
	    my $matchchain = $match[$i]{"chain"};
	    my $chainind = $chains{$matchchain};
	    my $outfile = $file;
	    $outfile =~ s/\.pdb/_align_$match[$i]{'ID'}_${matchchain}${chainind}\.pdb/;

	    open(MYOUT,">seqalign/$outfile");
	    open(MYIN,"$file");
	    while (my $line = <MYIN>) {
		if ((substr($line,0,4) eq "ATOM" || substr($line,0,4) eq "TER ") || substr($line,0,6) eq "HETATM") { 
		    my $chain = substr($line,21,1);
		    my $resnum = substr($line,22,4);
		    $resnum =~ s/\s+//;
		    if ($chain eq $matchchain) {  # write out modified residue number
			my $newnum;
			if ($resnum < $match[$i]{"startres"}) {
			    $newnum = $resnum - $match[$i]{"startres"} + $alignres;
			} else {  # shift up, taking dashes in sequence file into account
			    my $n = scalar @{$ndashbef{$match[$i]{'ID'}}};
			    my $tmp = $resnum-$match[$i]{"startres"};
			    if ($tmp >= $n) {
				$tmp = $n -1;
			    }
			    $newnum = $resnum - $match[$i]{"startres"} + $alignres + $ndashbef{$match[$i]{'ID'}}[$tmp];
			}
			if ($newnum < 0) {
			    if (-$newnum > $shiftupnexttime) {
				$shiftupnexttime = -$newnum;
			    }
			    $newnum = $resnum;
			    substr($line,21,1,"X");  # put residue instead on chain X, write a warning at the end
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

if ($shiftupnexttime > 0) {
    my $tmp = $shiftupnexttime + $alignres;
    print(STDERR "Warning:  negative residue indices were shifted to chain X. 
To avoid this next time, use alignres = $tmp\n");
}

sub getmatch {

# looks for a match between the two sequences
# returns an array of two integer arrays,integer reporting on the match strength, and the starting residue
#
#    strength 1 = matches, but with X residues
#    strength 2 = matches, but with missing residues
#    strength 3 = matches

    my $seqpdb = shift;
    my $seqalign = shift;

    my @matches = ();

    if (length($seqpdb) >= length($seqalign)) {
	my @pdbres = split(//,$seqpdb);
	my @alignres = split(//,$seqalign);
	for my $startres (0..length($seqpdb)-length($seqalign)) {
	    my $match = 3;
	    my $nresmatch = 0;
	    my $i = 0;
	    for my $res ($startres..$startres+length($seqalign)-1) {
		
		if ($pdbres[$res] eq "X") {
		    $match = 1;
		} elsif ($pdbres[$res] eq " ") {
		    $match = 2;
		} elsif ($pdbres[$res] ne $alignres[$i]) {
		    $match = 0;
		    last;
		} else {
		    $nresmatch++;
		}
		$i++;
	    }
	    my $fracmatch = $nresmatch/(length($seqalign));
	    if ($match > 0 && $fracmatch > 0.5) {
		my %tmp = ("strength" => $match, "startres" => $startres);
		push(@matches,\%tmp);
	    }
	}
    }

    return \@matches;		
}
