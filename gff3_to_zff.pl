#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use Getopt::Std;

die "usage: $0 <fasta file> <gff file>\n" unless @ARGV == 2;
my ($fasta, $gff) = @ARGV;

my %c;
my $gf;
if ($gff =~ /\.gz$/) {open($gf, "gunzip -c $gff |") or die}
else                 {open($gf, $gff) or die}

while (<$gf>) {
	next if /^#/;
	my ($chr, $m, $f, $beg, $end, $ph, $str, $fr, $grp) = split(/\t/, $_);
	next unless $f eq 'CDS';
	my ($name) = $grp =~ /Parent=(\S+)/;
	($beg, $end) = ($end, $beg) if $str eq '-';
	push(@{$c{$chr}{$name}}, [$beg, $end]);
}

my @id;
my $ff;
if ($fasta =~ /\.gz$/) {open($ff, "gunzip -c $fasta |") or die}
else                   {open($ff, $fasta) or die}
while (<$ff>) {
	chomp;
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


