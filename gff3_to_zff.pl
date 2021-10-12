#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';

my %c;
open(my $gf, "genes.gff") or die;
while (<$gf>) {	
	next if /^#/;
	my ($chr, $m, $f, $beg, $end, $ph, $str, $fr, $grp) = split(/\t/, $_);
	next unless $f eq 'CDS';
	my ($name) = $grp =~ /Parent=(\S+)/;
	($beg, $end) = ($end, $beg) if $str eq '-';
	push(@{$c{$chr}{$name}}, [$beg, $end]);
}

my @id;
open(my $ff, "gunzip -c genome.fa.gz |") or die;
while (<$ff>) {
	if (/>(\S+)/) {
		push @id, $1;
	}
}

foreach my $chr (@id) {
	print ">$chr\n";
	foreach my $gene (keys %{$c{$chr}}) {
		foreach my $exon (@{$c{$chr}{$gene}}) {
			print "Exon $exon->[0] $exon->[1] $gene\n";
		}
	}
}


