#!/usr/bin/perl
use strict; use warnings;

my %Feature = (
	'Exon'  => 'CDS',
	'Einit' => 'CDS',
	'Eterm' => 'CDS',
	'Esngl' => 'CDS',
);

print "##gff-version 3\n";
my $seq;
my %H;
while (<>) {
	if (/^>(\S+)/) {
		print "#region $1\n";
		$seq = $1;
	} else {
		my @f = split;
		if (@f == 4) {
			my $strand = $f[1] < $f[2] ? '+' : '-';
			if ($strand eq '-') {($f[1], $f[2]) = ($f[2], $f[1])}
			print join("\t", $seq, 'snap', $Feature{$f[0]}, $f[1], $f[2], '.',
				$strand, "Name=$f[3]"), "\n";
		} elsif (@f == 9) {
			print join("\t", $seq, 'snap', $Feature{$f[0]}, $f[1], $f[2], '.',
				$f[3], "Name=$f[8]"), "\n";
		} else {die "input does not appear to be ZFF"}
	}
}
__END__
