#!/usr/bin/perl
$| =1;
use strict;
use warnings;
use sigtrap;
use Getopt::Std;
use vars qw(
	$opt_r $opt_o $opt_x $opt_i $opt_e
	$opt_A $opt_D $opt_M $opt_S $opt_C $opt_I $opt_N
	$opt_a $opt_3 $opt_p $opt_5 $opt_t $opt_Z $opt_1 $opt_c
);
getopts('roxi:e:A:D:M:S:C:I:N:a:3:p5:tZ:1c:');

my $ACCEPTOR = "0:30";
my $DONOR    = "0:9";
my $START    = "0:15";
my $STOP     = "0:9";
my $CODING   = "4";
my $INTRON   = "4";
my $INTER    = "4";
my $PROM;
my $UTR5     = "";
my $UTR3     = "";
my $POLYA    = "";

my $UTR5Length; # defined below;
my $UTR3Length; # defined below
my $InterLength = 500;
my $EsnglLength = 1000;

die "
usage: hmm-assembler.pl <name> <directory of files from forge>
options:
  -i  <length>       [$InterLength]
  -e  <length>       [$EsnglLength]
  -A  <order:length> [$ACCEPTOR]
  -D  <order:length> [$DONOR]
  -M  <order:length> [$START]
  -S  <order:length> [$STOP]
  -C  <order>        [$CODING]
  -I  <order>        [$INTRON]
  -N  <order>        [$INTER]
  -3  <order:length> [$UTR3]  include 3'UTR model, requires -a
  -a  <order:length> [$POLYA]  include PolyA model, requires -3
  -5  <order:length> [$UTR5]  include 5'UTR moel, requires -p
  -p                     include generic promoter model, requires -5
  -r                     include generic repeat model
  -o                     include reverse ORF model
  -x                     use explicit duration intron model
  -t                     include C.elegans trans-splicing, requires -p, -5
  -Z  <clade>            sets clade-specific values (worm, fly, plant)
  -1                     single gene model
  -c <score>             include GC-AG splice donor model
" unless @ARGV == 2;


my ($NAME, $DIR) = @ARGV;

my $REPEATS      = $opt_r;
my $REVERSE_ORF  = $opt_o;
my $EXPLICIT     = $opt_x;
my $TRANS_SPLICE = $opt_t;
my $SPECIES      = $opt_Z? $opt_Z : "";
my $SINGLE_GENE  = $opt_1;
my $GC_AG        = $opt_c;
$ACCEPTOR = $opt_A if $opt_A;
$DONOR    = $opt_D if $opt_D;
$START    = $opt_M if $opt_M;
$STOP     = $opt_S if $opt_S;
$CODING   = $opt_C if $opt_C;
$INTRON   = $opt_I if $opt_I;
$INTER    = $opt_N if $opt_N;
$POLYA    = $opt_a if $opt_a;
$UTR3     = $opt_3 if $opt_3;
$PROM     = $opt_p;
$UTR5     = $opt_5 if $opt_5;
$InterLength = $opt_i if $opt_i;
$EsnglLength = $opt_e if $opt_e;

if (($POLYA and !$UTR3) or (!$POLYA and $UTR3)) {
	die "both -a and -3 must be specified";
}
if ($POLYA and $UTR3) {
	my ($order, $length) = split(/:/, $UTR3);
	$UTR3 = $order; # reassign
	$UTR3Length = $length; # geometric distribution
}

if (($PROM and !$UTR5) or (!$PROM and $UTR5)) {
	die "both -p and -5 must be specified";
}
if ($PROM and $UTR5) {
	my ($order, $length) = split(/:/, $UTR5);
	$UTR5 = $order;
	$UTR5Length = $length;
}

if ($TRANS_SPLICE) {
	die "both -p and -5 must be specified" unless $PROM and $UTR5;
}

####################
# Species override #
####################

if ($SPECIES eq 'worm') {
	$ACCEPTOR = '0:15';
	$REVERSE_ORF = 1;
	$GC_AG = -5;
} elsif ($SPECIES eq 'plant') {
	$ACCEPTOR = '0:20';
} elsif ($SPECIES eq 'fly') {
	$ACCEPTOR = '0:30';
	$REVERSE_ORF = 1;
} elsif ($SPECIES =~ /\S/) {
	die "unrecognized clade ($SPECIES)";
}

#####################
# Single gene model #
#####################
if ($SINGLE_GENE) {die "single gene not supported yet"}

##########
# States #
##########
my $States = 6;
my %State = (
	Einit  => {init => 0, term => 0, min => 3,   max => -1, dur => 'explicit'},
	Exon   => {init => 0, term => 0, min => 6,   max => -1, dur => 'explicit'},
	Eterm  => {init => 0, term => 0, min => 3,   max => -1, dur => 'explicit'},
	Esngl  => {init => 0, term => 0, min => 150, max => -1, dur => 'explicit'},
	Inter  => {init => 0.9, term => 0.9, min => 0, max => 0, dur => 'geometric'},
	Intron => {init => 0.1, term => 0.1, min => 0, max => 0, dur => 'geometric'},
);

if ($REPEATS) {
	$State{Repeat} = {init => 0, term => 0, min => 100, max => -1, dur => 'explicit'};
	$States++;
}

if ($REVERSE_ORF) {
	$State{ORF} = {init => 0, term => 0, min => 100, max => -1, dur => 'explicit'};
	$States++;
}

if ($POLYA) {
	$State{PolyA} = {init => 0,   term => 0,   min => 1, max => 1, dur => 'explicit'};
	$State{UTR3}  = {init => 0.1, term => 0.1, min => 0, max => 0, dur => 'geometric'};
	$States += 2;
}

if ($PROM) {
	$State{Prom} = {init => 0,   term => 0,   min => 1, max => 1, dur => 'explicit'};
	$State{UTR5} = {init => 0.1, term => 0.1, min => 0, max => 0, dur => 'geometric'};
	$States += 2;
}

if ($TRANS_SPLICE) {
	$State{TSS} = {init => 0, term => 0, min => 1, max => 1, dur => 'explicit'};
	$States += 1;
}

if ($EXPLICIT) {
	$State{Intron}{dur} = 'explicit';
	$State{Intron}{min} = 1;
	$State{Intron}{max} = -1;
}

###############
# Transitions #
###############
my $Transitions = 4;
my %Transition = (
	Einit => {Intron => 1},
	Esngl => {Inter => 1},
	Eterm => {Inter => 1},
	Exon  => {Intron => 1},
);

open(FILE, "$DIR/transitions");
while (<FILE>) {
	my ($s1, $s2, $prob) = split;
	$Transition{$s1}{$s2} = $prob;
	$Transitions++;
}
close FILE;

# optional section

if ($REPEATS) {
	$Transition{Intron}{Repeat} = 1;
	$Transition{Repeat}{Intron} = 1;
	$Transition{Inter}{Repeat} = 1;
	$Transition{Repeat}{Inter} = 1;
	$Transitions += 4;
}

if ($REVERSE_ORF) {
	$Transition{Intron}{ORF} = 1;
	$Transition{ORF}{Intron} = 1;
	$Transition{Inter}{ORF} = 1;
	$Transition{ORF}{Inter} = 1;
	$Transitions += 4;
}

if ($POLYA) {
	$Transition{Esngl} = {UTR3 => 1};
	$Transition{Eterm} = {UTR3 => 1};
	$Transition{UTR3}  = {PolyA => 1};
	$Transition{PolyA} = {Inter => 1};
	$Transitions += 2;
}

if ($PROM) {
	$Transition{Inter}{Prom} = 1;
	$Transition{Prom}{UTR5}  = 1;
	$Transition{UTR5}{Esngl} = $Transition{Inter}{Esngl};
	$Transition{UTR5}{Einit} = $Transition{Inter}{Einit};
	delete $Transition{Inter}{Esngl};
	delete $Transition{Inter}{Einit};
	$Transitions += 2;
}

if ($TRANS_SPLICE) {
	$Transition{Inter}{Prom} = 0.5;
	$Transition{Inter}{TSS}  = 0.5; # what are the real figures?
	$Transitions++;
}

###############
# Phase prefs #
###############
my $Phaseprefs = `cat $DIR/phaseprefs`;

#############
# Durations #
#############
my $Durations = 6;
my %Duration = (
	Einit => 1,
	Eterm => 1,
	Esngl => 1,
	Exon  => 1,
	Intron => 1,
	Inter => 1,
);
if ($REPEATS) {
	$Duration{Repeat} = 1;
	$Durations++;
}
if ($REVERSE_ORF) {
	$Duration{ORF} = 1;
	$Durations++;
}
if ($POLYA) {
	$Duration{PolyA} = 1;
	$Duration{UTR3}  = 1;
	$Durations += 2;
}
if ($PROM) {
	$Duration{Prom} = 1;
	$Duration{UTR5} = 1;
	$Durations += 2;
}
if ($TRANS_SPLICE) {
	$Duration{TSS} = 1;
	$Durations++;
}

##########
# Models #
##########
my $Models = 7;
my %Model = (
	Acceptor => $ACCEPTOR,
	Donor    => $DONOR,
	Start    => $START,
	Stop     => $STOP,
	Coding   => $CODING,
	Intron   => $INTRON,
	Inter    => $INTER,
);
if ($REPEATS) {
	$Model{Repeat} = 'Repeat';
	$Models++;
}
if ($POLYA) {
	$Model{PolyA} = $POLYA;
	$Model{UTR3}  = $UTR3;
	$Models += 2;
}
if ($PROM) {
	$Model{Prom} = 'Prom';
	$Model{UTR5} = $UTR5;
	$Models += 2;
}
if ($TRANS_SPLICE) {
	$Model{TSS} = 'TSS';
	$Models++;
}
# no need for ORF models, they are the same as coding

############
# Assemble #
############

# header
print "zoeHMM $NAME $States $Transitions $Durations $Models\n";

# states
print "\n<STATES>\n\n";
foreach my $name (sort keys %State) {
	print join("\t", $name, $State{$name}{init}, $State{$name}{term},
		$State{$name}{min}, $State{$name}{max}, $State{$name}{dur}), "\n";
}

# transitions
print "\n<STATE_TRANSITIONS>\n\n";
foreach my $s1 (sort keys %Transition) {
	foreach my $s2 (sort keys %{$Transition{$s1}}) {
		print join("\t", $s1, $s2, $Transition{$s1}{$s2}), "\n";
	}
}

# phaseprefs
print "\n<PHASE_PREFERENCES>\n\n";
print $Phaseprefs;

# durations
print "\n<STATE_DURATIONS>\n\n";
foreach my $name (sort keys %Duration) {
	if ($name eq 'Repeat') {
		print "Repeat 1\n\tCONSTANT 0 -1\n\t\t0\n";
	}
	elsif ($name eq 'ORF') {
		my $hack = `cat $DIR/Exon-explicit.duration`;
		$hack =~ s/Exon/ORF/;
		print $hack;
	}
	elsif ($name eq 'PolyA') {
		print "PolyA 1\n\tCONSTANT 0 -1\n\t\t0\n"; # may want to change that
	}
	elsif ($name eq 'Prom') {
		print "Prom 1\n\tCONSTANT 0 -1\n\t\t0\n"; # may want to change this too
	}
	elsif ($name eq 'TSS') {
		print "TSS 1\n\tCONSTANT 0 -1\n\t\t0\n"; # here as well
	}
	elsif ($name eq 'Inter') {
		print "Inter 1\n\tGEOMETRIC 0 -1\n\t\t$InterLength\n";
	}
	elsif ($name eq 'Esngl') {
		print "Esngl 1\n\tGEOMETRIC 0 -1\n\t\t$EsnglLength\n";
	}
	elsif ($name eq 'UTR3') {
		print "UTR3 1\n\tGEOMETRIC 0 -1\n\t\t$UTR3Length\n";
	}
	elsif ($name eq 'UTR5') {
		print "UTR5 1\n\tGEOMETRIC 0 -1\n\t\t$UTR5Length\n";
	}
	else {
		my $file = "$name-$State{$name}{dur}.duration";
		system("cat $DIR/$file");
	}
	print "\n";
}

# models

my $TSS_MODEL = 'TSS SDT 2 1 4 2 0.000
	AG WMM 15 11 4 0 0.000
		0.440	-0.800	-1.561	0.790	
		0.452	-0.664	-1.565	0.734	
		0.546	-0.660	-1.787	0.694	
		0.681	-1.267	-1.862	0.771	
		0.743	-1.716	-1.906	0.811	
		0.175	-1.667	-2.024	1.208	
		-2.278	-2.834	-4.043	1.845	
		-5.146	-4.228	-6.441	1.966	
		-1.458	-0.723	-1.433	1.411	
		-2.982	1.737	-7.579	-0.906	
		1.999	-9.900	-9.900	-9.900	
		-9.900	-9.900	1.999	-9.900	
		0.684	-0.753	0.351	-0.929	
		0.220	-0.506	-0.606	0.560	
		0.235	-0.166	-0.481	0.281	
	NN TRM 0 0 0 0 0.000
'; # C. elegans weight matrix for acceptor site

print "<SEQUENCE_MODELS>\n\n";
foreach my $name (sort keys %Model) {
	my ($order, $length);
	if ($Model{$name} eq 'Repeat') {
		print "Repeat LUT 1 0 5 0 0\n\t-1 -1 -1 -1 1\n\n";
		next;
	}
	elsif ($Model{$name} eq 'Prom') {
		print "Prom WMM 2 0 4 0 0\n\t0 0 0 0\n\t0 0 0 0\n\n";
		next;
	}
	elsif ($Model{$name} eq 'TSS') {
		print $TSS_MODEL, "\n";
		next;
	}
	elsif ($Model{$name} =~ /:/) {
		($order, $length) = split(/:/, $Model{$name});
	}
	else {
		if ($Model{$name} =~ /^(\d+)\+/) {
			$order = $1;
			$length = ($order +1) . "+";
		}
		else {
			$order = $Model{$name};
			$length = $order +1;
		}
	}
	my $file = "$DIR/$name-$order-$length.model";
	if ($name eq 'Donor' and $GC_AG) {
		print "Donor SDT 2 0 4 3 0.0\n";
		my @gc;
		open(DN, "Donor-0-9.model") or die;
		while (<DN>) {
			my @f = split;
			push @gc, \@f;
		}
		close DN;
		shift @gc; pop @gc; # unnecessary bits
		$gc[0][0] = 'GC';
		$gc[0][6] = $GC_AG; # assign score here
		($gc[5][1], $gc[5][3]) = ($gc[5][3], $gc[5][1]); # swap T and C values
		print "\t@{$gc[0]}\n";
		for (my $i = 1; $i < @gc; $i++) {
			print "\t\t", join("\t", @{$gc[$i]}), "\n";
		}
		open(IN, $file) or die;
		my $head = <IN>; # throw-away
		while (<IN>) {print}
	} else {
		system("cat $file");
	}
	print "\n";
}


# C. elegans trans-splice site
# using splice acceptor site weight matrix

__END__

Copyright (C) 2003-2004 Ian Korf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
