SNAP Documentation
==================

## Introduction ##

SNAP is a general purpose gene finding program suitable for both eukaryotic and
prokaryotic genomes. SNAP is an acroynm for Semi-HMM-based Nucleic Acid Parser.

## Reference ##

Korf I. Gene finding in novel Genomes. BMC Bioinformatics 2004, 5:59

## Contact ##

I appreciate bug reports, comments, and suggestions. My current contact
information is:

* email: ifkorf@ucdavis.edu
* web: http://korflab.ucdavis.edu

## License ##

This software is covered by the MIT License.

## Files and Directories ##

    DNA               Contains some sample sequences
    HMM               Contains SNAP parameter files
    LICENSE           MIT license
    Makefile          For compiling
    Makefile.include  Automatically generated, should not be edited
    Zoe               A library containing lots of the base functions
    fathom.c          Utility for investigating sequences and annotation
    forge.c           Parameter estimation
    hmm-assembler.pl  Creates HMMs for SNAP
    snap.c            Gene prediction program

## Your favorite genome... ##

If you wish to train SNAP for a new genome, please contact me. Parameter
estimation is not particularly difficult, but the procedure is not well
documented and I have only included the minimum applications here. I've included
the basic strategy at the end of this document.


INSTALLATION INSTRUCTIONS
=========================

The software is routinely compiled and tested on Mac OS X. It should compile
fine on any Linux/Unix type operating systems but new compilers sometimes
complain, so please let me know if you have problems.

## Enviroment ##

The ZOE environment variable can be  used by SNAP to find the HMM files. Set
this to the directory containing this file. For example, if you unpackaged the
tar-ball in /usr/local, set the ZOE environment variable to /usr/local/Zoe
    
	setenv ZOE /usr/local/Zoe # csh, tcsh, etc
	export ZOE=/usr/local/Zoe # sh, bash, etc

If you do not use the ZOE environment variable, you can still use
SNAP but you must specify the explict path to the parameter file.

## Compiling ##

Provided you have gcc and make, compiling should be as simple as:

	make

Testing

	./snap HMM/thale DNA/thale.dna.gz
	./snap HMM/worm DNA/worm.dna.gz


PARAMETER ESTIMATION
====================

Sequences must be in FASTA format. It's a good idea if you don't have genes that
are too related to each other.

Gene structures must be in ZFF format. What is ZFF? It is a non-standard format
(ie. nobody uses it but me) that bears resemblence to FASTA and GFF (both true
standards). There are two styles of ZFF, the short format and the long format.
In both cases, the sequence records are separated by a definition line, just
like FASTA. In the short format, there are 4 fields: Label, Begin, End, Group.
The 4th field is optional. Label is a controlled vocabulary (see zoeFeature.h
for a complete list). All exons of a gene (or more appropriately a
transcriptional unit) must share the same unique group name. The strand of the
feature is implied in the coordinates, so if Begin > End, the feature is on the
minus strand. Here's and example of the short format with two sequences, each
containing a single gene on the plus strand:

    >sequence-1
    Einit    201    325   Y73E7A.6
    Eterm   2175   2319   Y73E7A.6
    >sequence-2
    Einit    201    462   Y73E7A.7
    Exon    1803   2031   Y73E7A.7
    Exon    2929   3031   Y73E7A.7
    Exon    3467   3624   Y73E7A.7
    Exon    4185   4406   Y73E7A.7
    Eterm   5103   5280   Y73E7A.7

The long format adds 5 fields between the coordinates and the group: Strand,
Score, 5'-overhang, 3'-overhang, and Frame. Strand is +/-. Score is any floating
point value. 5'- and 3'-overhang are the number of bp of an incomplete codon at
each end of an exon. Frame is the reading frame (0..2 and *not* 1..3). Here's an
example of the long format:

    >Y73E7A.6
    Einit    201    325   +    90   0   2   1   Y73E7A.6
    Eterm   2175   2319   +   295   1   0   2   Y73E7A.6
    >Y73E7A.7
    Einit    201    462   +   263   0   1   1   Y73E7A.7
    Exon    1803   2031   +   379   2   2   0   Y73E7A.7
    Exon    2929   3031   +   236   1   0   0   Y73E7A.7
    Exon    3467   3624   +   152   0   2   0   Y73E7A.7
    Exon    4185   4406   +   225   1   2   2   Y73E7A.7
    Eterm   5103   5280   +    46   1   0   2   Y73E7A.7

The most important part of parameter estimation is preparing a training set.
There are many ways to go about this. At the end, you want these in the ZFF
short format. Save the ZFF as genome.ann and the FASTA as genome.dna. The first
step is to look at some features of the genes:

    fathom genome.ann genome.dna -gene-stats 

Next, you want to verify that the genes have no obvious errors:

    fathom genome.ann genome.dna -validate

You may find some errors and warnings. Check these out in some kind of genome
browser and remove those that are real errors. Next, break up the sequences into
fragments with one gene per sequence with the following command:

    fathom genome.ann genome.dna -categorize 1000

There will be up to 1000 bp on either side of the genes. You will find
several new files.

    alt.ann, alt.dna (genes with alternative splicing)
    err.ann, err.dna (genes that have errors)
    olp.ann, olp.dna (genes that overlap other genes)
    wrn.ann, wrn.dna (genes with warnings)
    uni.ann, uni.dna (single gene per sequence)

Convert the uni genes to plus stranded with the command:

    fathom uni.ann uni.dna -export 1000 -plus

You will find 4 new files:

    export.aa   proteins corresponding to each gene
    export.ann  gene structure on the plus strand
    export.dna  DNA of the plus strand
    export.tx   transcripts for each gene

The parameter estimation program, forge, creates a lot of files. You probably
want to create a directory to keep things tidy before you execute the program.

    mkdir params
    cd params
    forge ../export.ann ../export.dna
    cd ..

Last is to build an HMM.

    hmm-assembler.pl my-genome params > my-genome.hmm

There are a number of options for forge and hmm-assembler.pl that I have not
described here. Hopefully I'll document these one day.
