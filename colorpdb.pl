#!/usr/bin/env perl

use strict;

use lib "/Users/alexrd/perl";
use mj qw ( %names %symb );

my $pdbfile = shift(@ARGV);
my $contactfile = shift(@ARGV);
my $seqfile = shift(@ARGV);

if (!defined $contactfile) {
    $contactfile = "contact_freq.dat";
}

if (!defined $seqfile) {
    $seqfile = "clustal_frag_align.dat";
}


sub help {
    print("
colorpdb.pl

Usage:  ./colorpdb.pl pdb_file contact_table (default: contact_freq.dat) sequence_align_file (default: clustal_frag_align.dat)

Alex Dickson
University of Michigan

This script uses a specified pdb as a template, and creates a set of pdbs with
contact frequency information in the beta factor column.  Given NAME.pdb, it will output:

NAME_prot.pdb
NAME_phob.pdb
NAME_phil.pdb
NAME_water.pdb

which show the frequency of contacts with all protein amino acids, hydrophobic amino
acids, hydrophilic amino acids, and water, respectively.

This also prints out pdbs that are colored by sequence conservation:

NAME_cons_exact.pdb  -- printing the frequency of exact matches
NAME_cons_type.pdb   -- printing the frequency of matches by type (hydrophobic, aromatic, +ve, -ve, etc.)

It uses the both the NAME.pdb file in the current directory, as well as the aligned pdb
file in the seqalign/ directory to get the amino acid numbering after alignment.
");
    exit(0);
}

my $flag = shift(@ARGV);
if ($flag eq "-h" || $flag eq "-help") {
    &help;
}

if (!-f $pdbfile) {    
    print("Error:  pdbfile does not exist!\n");
    &help;
}

(my $pdb) = ($pdbfile =~ m/([0-9A-Za-z]+)\.pdb/);
my $tmp = `ls seqalign/${pdb}*`;
my @alignfiles = split(/\n/,$tmp);

if (scalar @alignfiles == 0) {
    print("Error:  no aligned files in the seqalign/ directory!\n");
    &help;
}

my %alignseq;
# read alignfile
if (!-f $seqfile) {
    print("Error!  seqfile $seqfile is not defined!\n");
    &help;
}
open(MYIN,$seqfile);
while (my $line = <MYIN>) {
    chomp($line);
    my @arr = split(/\s+/,$line);
    if (defined $arr[0] && $arr[0] ne "CLUSTAL") {
	$alignseq{$arr[0]} .= "$arr[1]";
    }
}
close(MYIN);

(my $exactref, my $typeref) = &getcons(\%alignseq);
my @match_exact = @{$exactref};
my @match_type = @{$typeref};

my @labels;
my %resinfo;
my $firstres = undef;
# read in the table
open(MYIN,"$contactfile");
my $line = <MYIN>;
chomp($line);
my @labels = split(/ /,$line);
while ($line = <MYIN>) {
    chomp($line);
    my @arr = split(/ /,$line);
    if (!defined $firstres) {
	$firstres = $arr[0];
    }
    for my $i (1..@arr-1) {
	$resinfo{$arr[0]}{$labels[$i]} = $arr[$i];
    }
}

my @tk = keys %alignseq;
my $l = length($alignseq{$tk[0]});
my $lastres = $firstres + $l - 1;

# get avg properties:  prot, phob, phil, water
my @phob = qw (val leu ile phe met);
my @phil = qw (asp asn glu gln ser lys arg his thr);
my @others = qw (gly ala cys tyr trp pro);
for my $res (keys %resinfo) {
    $resinfo{$res}{"phob"} = 0;
    $resinfo{$res}{"prot"} = 0;
    $resinfo{$res}{"phil"} = 0;
    $resinfo{$res}{"water"} = 0;
    for my $lab (keys %{$resinfo{$res}}) {
	for (@phob) {
	    if ($lab eq $_) {
		$resinfo{$res}{"phob"} += $resinfo{$res}{$_};
		$resinfo{$res}{"prot"} += $resinfo{$res}{$_};
	    }
	}

	for (@phil) {
	    if ($lab eq $_) {
		$resinfo{$res}{"phil"} += $resinfo{$res}{$_};
		$resinfo{$res}{"prot"} += $resinfo{$res}{$_};
	    }
	}

	for (@others) {
	    if ($lab eq $_) {
		$resinfo{$res}{"prot"} += $resinfo{$res}{$_};
	    }
	}

	if ($lab eq "hoh") {
	    $resinfo{$res}{"water"} += $resinfo{$res}{$_};
	}
    }
}

# read through pdb, printing out new ones as you go
if (!-d "color") {
    mkdir("color");
}
for my $i (0..@alignfiles-1) {
    (my $alignchain) = ($alignfiles[$i] =~ m/([A-Z])\.pdb/);
    open(MYIN,"$alignfiles[$i]");
    my $out = $alignfiles[$i];
    $out =~ s/seqalign/color/;
    my $phobout = $out;
    $phobout =~ s/align/phob/;
    my $philout = $out;
    $philout =~ s/align/phil/;
    my $protout = $out;
    $protout =~ s/align/prot/;
    my $waterout = $out;
    $waterout =~ s/align/water/;
    my $exactout = $out;
    $exactout =~ s/align/match_exact/;
    my $typeout = $out;
    $typeout =~ s/align/match_type/;
    open(PHOBOUT,">$phobout");
    open(PHILOUT,">$philout");
    open(PROTOUT,">$protout");
    open(WATEROUT,">$waterout");
    open(EXACTOUT,">$exactout");
    open(TYPEOUT,">$typeout");

    while (my $line = <MYIN>) {
	if (substr($line,0,4) eq "ATOM" || substr($line,0,6) eq "HETATM") {
	    my $chain = substr($line,21,1);
	    my $resnum = substr($line,22,4);
	    $resnum =~ s/\s+//;
	    
	    if ($chain eq "$alignchain" && exists $resinfo{$resnum}) {
		my $tline = $line;
		substr($tline,55,5,sprintf("%5.2f",$resinfo{$resnum}{"phob"}));
		print(PHOBOUT "$tline");

		$tline = $line;
		substr($tline,55,5,sprintf("%5.2f",$resinfo{$resnum}{"phil"}));
		print(PHILOUT "$tline");

		$tline = $line;
		substr($tline,55,5,sprintf("%5.2f",$resinfo{$resnum}{"prot"}));
		print(PROTOUT "$tline");

		$tline = $line;
		substr($tline,55,5,sprintf("%5.2f",$resinfo{$resnum}{"water"}));
		print(WATEROUT "$tline");

		$tline = $line;		
		substr($tline,55,5,sprintf("%5.2f",$match_exact[$resnum-$firstres]));
		print(EXACTOUT "$tline");

		$tline = $line;		
		substr($tline,55,5,sprintf("%5.2f",$match_type[$resnum-$firstres]));
		print(TYPEOUT "$tline");

	    } else {
		substr($line,55,5,sprintf("%5.2f",0));
		print(PHOBOUT "$line");
		print(PHILOUT "$line");
		print(PROTOUT "$line");
		print(WATEROUT "$line");
		print(EXACTOUT "$line");
		print(TYPEOUT "$line");
	    }
	}
    }
    close(MYIN);
    close(PHOBOUT);
    close(PHILOUT);
    close(PROTOUT);
    close(WATEROUT);
    close(EXACTOUT);
    close(TYPEOUT);

# visualize first file
    my $vmd = "/Developer/Applications/VMD\\ 1.9.1.app/Contents/MacOS/startup.command";
    system("$vmd -e viscolors.tcl -args $phobout $philout $protout $waterout $exactout $typeout $firstres $lastres $alignchain");
}




sub getcons {
    my %seq = %{@_[0]};
    my @keys = keys %seq;

    my $n = length($seq{$keys[0]});

    my @exact;
    my @type;

    for my $i (0..$n-1) {
	my %aas;
	my %aas_type;

	my $nused = 0;

	for my $key (@keys) {
	    my $aa = substr($seq{$key},$i,1);
	    if ($aa ne "-") {
		my $type = &gettype($aa);
		if (!exists $aas{$aa}) {
		    $aas{$aa} = 0;
		}
		$aas{$aa}++;
		if (!exists $aas_type{$type}) {
		    $aas_type{$type} = 0;
		}
		$aas_type{$type}++;
		$nused++;
	    }
	}

	my $norm = $nused*($nused-1);

	if ($norm > 0) {
	    my $sum = 0;
	    for (keys %aas) {
		$sum += ($aas{$_}-1)*$aas{$_};
	    }
	    $exact[$i] = $sum/$norm;
	    
	    $sum = 0;
	    for (keys %aas_type) {
		$sum += ($aas_type{$_}-1)*$aas_type{$_};
	    }
	    $type[$i] = $sum/$norm;
	} else {
	    $exact[$i] = 1;
	    $type[$i] = 1;
	}
    }

    return(\@exact,\@type);
}

sub gettype {

# uses optimal amino acid grouping predicted by maximum variance (see Table I of Wrabl and Grishin, Proteins, (61) 523-534, 2005)

    my $aa = lc(shift);
    if ("wfy" =~ m/$aa/) {
	return 0;
    } elsif ("mliv" =~ m/$aa/) {
	return 1;
    } elsif ("c" =~ m/$aa/) {
	return 2;
    } elsif ("g" =~ m/$aa/) {
	return 3;
    } elsif ("p" =~ m/$aa/) {
	return 4;
    } elsif ("ats" =~ m/$aa/) {
	return 5;
    } elsif ("nde" =~ m/$aa/) {
	return 6;
    } elsif ("hqrk" =~ m/$aa/) {
	return 7;
    } else {
	die("Error! amino acid $aa has no type");
    }

    
    
}
    
	
