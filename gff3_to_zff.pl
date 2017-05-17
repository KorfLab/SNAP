#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';

my %c;
while (<>) {
	my ($chr, $m, $f, $beg, $end, $ph, $str, $fr, $grp) = split;
	next unless $f eq 'CDS';
	my ($name) = $grp =~ /Parent=(\S+)/;
}
