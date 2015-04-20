#!/usr/bin/env perl

package mj;
use strict;
use Exporter;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(%mj %names %symb);
%EXPORT_TAGS = ( DEFAULT => [qw(%mj %names %symb)] );

our %mj;

$mj{"cys"}{"cys"} = -5.44;

$mj{"met"}{"cys"} = -4.99;
$mj{"met"}{"met"} = -5.46;

$mj{"phe"}{"cys"} = -5.80;
$mj{"phe"}{"met"} = -6.56;
$mj{"phe"}{"phe"} = -7.26;

$mj{"ile"}{"cys"} = -5.50;
$mj{"ile"}{"met"} = -6.02;
$mj{"ile"}{"phe"} = -6.84;
$mj{"ile"}{"ile"} = -6.54;

$mj{"leu"}{"cys"} = -5.83;
$mj{"leu"}{"met"} = -6.41;
$mj{"leu"}{"phe"} = -7.28;
$mj{"leu"}{"ile"} = -7.04;
$mj{"leu"}{"leu"} = -7.37;

$mj{"val"}{"cys"} = -4.96;
$mj{"val"}{"met"} = -5.32;
$mj{"val"}{"phe"} = -6.29;
$mj{"val"}{"ile"} = -6.05;
$mj{"val"}{"leu"} = -6.48;
$mj{"val"}{"val"} = -5.52;

$mj{"trp"}{"cys"} = -4.95;
$mj{"trp"}{"met"} = -5.55;
$mj{"trp"}{"phe"} = -6.16;
$mj{"trp"}{"ile"} = -5.78;
$mj{"trp"}{"leu"} = -6.14;
$mj{"trp"}{"val"} = -5.18;
$mj{"trp"}{"trp"} = -5.06;

$mj{"tyr"}{"cys"} = -4.16;
$mj{"tyr"}{"met"} = -4.91;
$mj{"tyr"}{"phe"} = -5.66;
$mj{"tyr"}{"ile"} = -5.25;
$mj{"tyr"}{"leu"} = -5.67;
$mj{"tyr"}{"val"} = -4.62;
$mj{"tyr"}{"trp"} = -4.66;
$mj{"tyr"}{"tyr"} = -4.17;

$mj{"ala"}{"cys"} = -3.57;
$mj{"ala"}{"met"} = -3.94;
$mj{"ala"}{"phe"} = -4.81;
$mj{"ala"}{"ile"} = -4.58;
$mj{"ala"}{"leu"} = -4.91;
$mj{"ala"}{"val"} = -4.04;
$mj{"ala"}{"trp"} = -3.82;
$mj{"ala"}{"tyr"} = -3.36;
$mj{"ala"}{"ala"} = -2.72;

$mj{"gly"}{"cys"} = -3.16;
$mj{"gly"}{"met"} = -3.39;
$mj{"gly"}{"phe"} = -4.13;
$mj{"gly"}{"ile"} = -3.78;
$mj{"gly"}{"leu"} = -4.16;
$mj{"gly"}{"val"} = -3.38;
$mj{"gly"}{"trp"} = -3.42;
$mj{"gly"}{"tyr"} = -3.01;
$mj{"gly"}{"ala"} = -2.31;
$mj{"gly"}{"gly"} = -2.24;

$mj{"thr"}{"cys"} = -3.11;
$mj{"thr"}{"met"} = -3.51;
$mj{"thr"}{"phe"} = -4.28;
$mj{"thr"}{"ile"} = -4.03;
$mj{"thr"}{"leu"} = -4.34;
$mj{"thr"}{"val"} = -3.46;
$mj{"thr"}{"trp"} = -3.22;
$mj{"thr"}{"tyr"} = -3.01;
$mj{"thr"}{"ala"} = -2.32;
$mj{"thr"}{"gly"} = -2.08;
$mj{"thr"}{"thr"} = -2.12;

$mj{"ser"}{"cys"} = -2.86;
$mj{"ser"}{"met"} = -3.03;
$mj{"ser"}{"phe"} = -4.02;
$mj{"ser"}{"ile"} = -3.52;
$mj{"ser"}{"leu"} = -3.92;
$mj{"ser"}{"val"} = -3.05;
$mj{"ser"}{"trp"} = -2.99;
$mj{"ser"}{"tyr"} = -2.78;
$mj{"ser"}{"ala"} = -2.01;
$mj{"ser"}{"gly"} = -1.82;
$mj{"ser"}{"thr"} = -1.96;
$mj{"ser"}{"ser"} = -1.67;

$mj{"asn"}{"cys"} = -2.59;
$mj{"asn"}{"met"} = -2.95;
$mj{"asn"}{"phe"} = -3.75;
$mj{"asn"}{"ile"} = -3.24;
$mj{"asn"}{"leu"} = -3.74;
$mj{"asn"}{"val"} = -2.83;
$mj{"asn"}{"trp"} = -3.07;
$mj{"asn"}{"tyr"} = -2.76;
$mj{"asn"}{"ala"} = -1.84;
$mj{"asn"}{"gly"} = -1.74;
$mj{"asn"}{"thr"} = -1.88;
$mj{"asn"}{"ser"} = -1.58;
$mj{"asn"}{"asn"} = -1.68;

$mj{"gln"}{"cys"} = -2.85;
$mj{"gln"}{"met"} = -3.30;
$mj{"gln"}{"phe"} = -4.10;
$mj{"gln"}{"ile"} = -3.67;
$mj{"gln"}{"leu"} = -4.04;
$mj{"gln"}{"val"} = -3.07;
$mj{"gln"}{"trp"} = -3.11;
$mj{"gln"}{"tyr"} = -2.97;
$mj{"gln"}{"ala"} = -1.89;
$mj{"gln"}{"gly"} = -1.66;
$mj{"gln"}{"thr"} = -1.90;
$mj{"gln"}{"ser"} = -1.49;
$mj{"gln"}{"asn"} = -1.71;
$mj{"gln"}{"gln"} = -1.54;

$mj{"asp"}{"cys"} = -2.41;
$mj{"asp"}{"met"} = -2.57;
$mj{"asp"}{"phe"} = -3.48;
$mj{"asp"}{"ile"} = -3.17;
$mj{"asp"}{"leu"} = -3.40;
$mj{"asp"}{"val"} = -2.48;
$mj{"asp"}{"trp"} = -2.84;
$mj{"asp"}{"tyr"} = -2.76;
$mj{"asp"}{"ala"} = -1.70;
$mj{"asp"}{"gly"} = -1.59;
$mj{"asp"}{"thr"} = -1.80;
$mj{"asp"}{"ser"} = -1.63;
$mj{"asp"}{"asn"} = -1.68;
$mj{"asp"}{"gln"} = -1.46;
$mj{"asp"}{"asp"} = -1.21;

$mj{"glu"}{"cys"} = -2.27;
$mj{"glu"}{"met"} = -2.89;
$mj{"glu"}{"phe"} = -3.56;
$mj{"glu"}{"ile"} = -3.27;
$mj{"glu"}{"leu"} = -3.59;
$mj{"glu"}{"val"} = -2.67;
$mj{"glu"}{"trp"} = -2.99;
$mj{"glu"}{"tyr"} = -2.79;
$mj{"glu"}{"ala"} = -1.51;
$mj{"glu"}{"gly"} = -1.22;
$mj{"glu"}{"thr"} = -1.74;
$mj{"glu"}{"ser"} = -1.48;
$mj{"glu"}{"asn"} = -1.51;
$mj{"glu"}{"gln"} = -1.42;
$mj{"glu"}{"asp"} = -1.02;
$mj{"glu"}{"glu"} = -0.91;

$mj{"his"}{"cys"} = -3.60;
$mj{"his"}{"met"} = -3.98;
$mj{"his"}{"phe"} = -4.77;
$mj{"his"}{"ile"} = -4.14;
$mj{"his"}{"leu"} = -4.54;
$mj{"his"}{"val"} = -3.58;
$mj{"his"}{"trp"} = -3.98;
$mj{"his"}{"tyr"} = -3.52;
$mj{"his"}{"ala"} = -2.41;
$mj{"his"}{"gly"} = -2.15;
$mj{"his"}{"thr"} = -2.42;
$mj{"his"}{"ser"} = -2.11;
$mj{"his"}{"asn"} = -2.08;
$mj{"his"}{"gln"} = -1.98;
$mj{"his"}{"asp"} = -2.32;
$mj{"his"}{"glu"} = -2.15;
$mj{"his"}{"his"} = -3.05;

$mj{"arg"}{"cys"} = -2.57;
$mj{"arg"}{"met"} = -3.12;
$mj{"arg"}{"phe"} = -3.98;
$mj{"arg"}{"ile"} = -3.63;
$mj{"arg"}{"leu"} = -4.03;
$mj{"arg"}{"val"} = -3.07;
$mj{"arg"}{"trp"} = -3.41;
$mj{"arg"}{"tyr"} = -3.16;
$mj{"arg"}{"ala"} = -1.83;
$mj{"arg"}{"gly"} = -1.72;
$mj{"arg"}{"thr"} = -1.90;
$mj{"arg"}{"ser"} = -1.62;
$mj{"arg"}{"asn"} = -1.64;
$mj{"arg"}{"gln"} = -1.80;
$mj{"arg"}{"asp"} = -2.29;
$mj{"arg"}{"glu"} = -2.27;
$mj{"arg"}{"his"} = -2.16;
$mj{"arg"}{"arg"} = -1.55;

$mj{"lys"}{"cys"} = -1.95;
$mj{"lys"}{"met"} = -2.48;
$mj{"lys"}{"phe"} = -3.36;
$mj{"lys"}{"ile"} = -3.01;
$mj{"lys"}{"leu"} = -3.37;
$mj{"lys"}{"val"} = -2.49;
$mj{"lys"}{"trp"} = -2.69;
$mj{"lys"}{"tyr"} = -2.60;
$mj{"lys"}{"ala"} = -1.31;
$mj{"lys"}{"gly"} = -1.15;
$mj{"lys"}{"thr"} = -1.31;
$mj{"lys"}{"ser"} = -1.05;
$mj{"lys"}{"asn"} = -1.21;
$mj{"lys"}{"gln"} = -1.29;
$mj{"lys"}{"asp"} = -1.68;
$mj{"lys"}{"glu"} = -1.80;
$mj{"lys"}{"his"} = -1.35;
$mj{"lys"}{"arg"} = -0.59;
$mj{"lys"}{"lys"} = -0.12;

$mj{"pro"}{"cys"} = -3.07;
$mj{"pro"}{"met"} = -3.45;
$mj{"pro"}{"phe"} = -4.25;
$mj{"pro"}{"ile"} = -3.76;
$mj{"pro"}{"leu"} = -4.20;
$mj{"pro"}{"val"} = -3.32;
$mj{"pro"}{"trp"} = -3.73;
$mj{"pro"}{"tyr"} = -3.19;
$mj{"pro"}{"ala"} = -2.03;
$mj{"pro"}{"gly"} = -1.87;
$mj{"pro"}{"thr"} = -1.90;
$mj{"pro"}{"ser"} = -1.57;
$mj{"pro"}{"asn"} = -1.53;
$mj{"pro"}{"gln"} = -1.73;
$mj{"pro"}{"asp"} = -1.33;
$mj{"pro"}{"glu"} = -1.26;
$mj{"pro"}{"his"} = -2.25;
$mj{"pro"}{"arg"} = -1.70;
$mj{"pro"}{"lys"} = -0.97;
$mj{"pro"}{"pro"} = -1.75;

our %names;
$names{"A"} = "ala";
$names{"R"} = "arg";
$names{"N"} = "asn";
$names{"D"} = "asp";
$names{"C"} = "cys";
$names{"Q"} = "gln";
$names{"E"} = "glu";
$names{"G"} = "gly";
$names{"H"} = "his";
$names{"I"} = "ile";
$names{"L"} = "leu";
$names{"K"} = "lys";
$names{"M"} = "met";
$names{"F"} = "phe";
$names{"P"} = "pro";
$names{"S"} = "ser";
$names{"T"} = "thr";
$names{"W"} = "trp";
$names{"Y"} = "tyr";
$names{"V"} = "val";

our %symb;
$symb{"ala"} = "A";
$symb{"arg"} = "R";
$symb{"asn"} = "N";
$symb{"asp"} = "D";
$symb{"cys"} = "C";
$symb{"gln"} = "Q";
$symb{"glu"} = "E";
$symb{"gly"} = "G";
$symb{"his"} = "H";
$symb{"ile"} = "I";
$symb{"leu"} = "L";
$symb{"lys"} = "K";
$symb{"met"} = "M";
$symb{"mse"} = "M";  # modified MET
$symb{"phe"} = "F";
$symb{"pro"} = "P";
$symb{"ser"} = "S";
$symb{"thr"} = "T";
$symb{"trp"} = "W";
$symb{"tyr"} = "Y";
$symb{"yof"} = "Y";  # modified TYR
$symb{"val"} = "V";

1;
