#!/usr/bin/env perl

use strict;
use File::Copy;

my @index = split(/:/,shift(@ARGV));
my $vmdtext = "/Developer/Applications/VMD\\ 1.9.1.app/Contents/MacOS/startup.command -dispdev text";

sub help {
    print("
alignstructures.pl

Alex Dickson
University of Michigan

This script takes a folder full of pdbs, and aligns them in space using VMD.
It assumes that alignbypattern.pl was used to align the pdbs by sequence, and 
looks for pdb files of the form NAME_align_X1.pdb, where X is the chain that is aligned.

The resulting structure-aligned pdbs are saved as NAME_struct_align_X1.pdb (etc..),
and the RMSD for each is written to rmsd_struct_align.dat.

The current vmd command (in text mode) is set to:
$vmdtext
");
    exit(0);
}

if (scalar @index == 0) {    	  
    print("Usage:  ind1:ind2:ind3:..\n");
    &help;
}

my @files = split(/\n/,`ls *align*.pdb`);

if (scalar @files == 0) {
    print("No aligned pdbs in this directory!\n");
    &help;
}

my $template = shift(@files);
(my $tmplchain) = ($template =~ m/align_([A-Z])/);

# write selection strings (with dashes in place of spaces, so its treated as one word)
my $sel = "resid-";
for (@index) { $sel .= "${_}-" ; } 
my $sel2 = $sel."and-backbone-and-chain-$tmplchain";

open(RMSDOUT,">rmsd_struct_align.dat");
print(RMSDOUT "$template 0\n");

for my $file (@files) {

    # amend selection string with the proper chain
    (my $chain) = ($file =~ m/align_([A-Z])/);
    my $sel1 .= $sel."and-backbone-and-chain-$chain";

    # align file using VMD
    system("$vmdtext -eofexit -args $file $template $sel1 $sel2 trmsd.out < align_2structs.tcl > t.log");

    # move t.pdb file to its unique name
    my $newfile = $file;
    $newfile =~ s/align/struct_align/;
    move("t.pdb",$newfile);

    # print rmsd
    my $rmsd = `cat trmsd.out`;
    print(RMSDOUT "$file $rmsd");

}
close(RMSDOUT);
