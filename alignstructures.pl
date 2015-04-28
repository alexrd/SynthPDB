#!/usr/bin/env perl

use strict;
use File::Copy;

my $vmdtext = "/Developer/Applications/VMD\\ 1.9.1.app/Contents/MacOS/startup.command -dispdev text";

my $alignfile = shift(@ARGV);
my $alignres = shift(@ARGV);
my $template = shift(@ARGV);

sub help {
    print("
alignstructures.pl

Usage:  alignstructures.pl alignfile alignres template.pdb

Alex Dickson
University of Michigan

This script takes a folder full of pdbs, and aligns them in space using VMD.
It assumes that groupalign.pl was used to align the pdbs by sequence, and 
looks for pdb files in the seqalign/ directory, which have the form *_Xn.pdb
where X is the chain that is aligned, and n is an index that keeps track of how many such pdbs 
there are.

The resulting structure-aligned pdbs are saved in the structalign/ directory as NAME_struct_align_X1.pdb (etc..),
and the RMSD for each is written to structalign/rmsd_struct_align.dat.

The first argument is the alignment file used in groupalign.pl

The second argument is the index of alignment used in groupalign.pl.

The third argument is the name of the pdb to be used as an alignment template, this should be one of the
files in the seqalign/ directory.

The current vmd command (in text mode) is set to:
$vmdtext
");
    exit(0);
}

if (!-f $alignfile) {    
    print("Error:  alignfile does not exist!\n");
    &help;
}

my @files = split(/\n/,`ls seqalign/*.pdb`);

if (scalar @files == 0) {
    print("Error:  no pdbs in the seqalign/ directory!\n");
    &help;
}

if (!-f $template) {
    print("Error:  template PDB does not exist!\n");
    &help;
}

(my $tmplchain) = ($template =~ m/([A-Z])[0-9]\.pdb/);
(my $tmplID) = ($template =~ m/align_([A-Z0-9]+_[0-9]+)_/);

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

if (!-d "structalign") {
    mkdir("structalign");
}

open(RMSDOUT,">structalign/rmsd_struct_align.dat");

for my $file (@files) {

    # get ID and chain of file
    (my $fileID) = ($file =~ m/align_([A-Z0-9]+_[0-9]+)_/);
    (my $chain) = ($file =~ m/([A-Z])[0-9]\.pdb/);

    # get list of common residues using alignseq
    my @common;
    for my $i (0..length($alignseq{$fileID})-1) {
	if (substr($alignseq{$fileID},$i,1) ne "-" && substr($alignseq{$tmplID},$i,1) ne "-") {
	    push(@common,$i+$alignres);
	}
    }

    if (scalar @common == 0) {
	print(STDERR "Error! no common residues between $fileID and $tmplID!");
    } else {
    
	# write selection strings (with dashes in place of spaces, so its treated as one word)
	my $sel = "resid-";
	for (@common) { $sel .= "${_}-" ; } 
	my $sel2 = $sel."and-name-CA-and-chain-$tmplchain";
	my $sel1 = $sel."and-name-CA-and-chain-$chain";
	
	# align file using VMD
	system("$vmdtext -eofexit -args $file $template $sel1 $sel2 trmsd.out < align_2structs.tcl > t.log");

	# move t.pdb file to its unique name
	my $newfile = $file;
	$newfile =~ s/seqalign/structalign/;
	$newfile =~ s/_align_/_struct_align_/;
	move("t.pdb",$newfile);

	# print rmsd
	my $rmsd = `cat trmsd.out`;
	if ($rmsd =~ m/X/) {
	    print(STDERR "Error in rmsd calculation between $fileID and $tmplID");
	} else {
	    print(RMSDOUT "$file $rmsd");
	}
    }
}
close(RMSDOUT);
