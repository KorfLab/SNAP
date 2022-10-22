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
documented and I have only included the minimum applications here. I've 
included the basic strategy at the end of this document.


INSTALLATION INSTRUCTIONS
=========================

The software is routinely compiled and tested on Mac and Linux on a variety of 
architectures. It should compile cleanly on any Linux/Unix type operating 
systems but new compilers sometimes complain, so please let me know if you have 
problems.

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

Sequences must be in FASTA format. It's a good idea if you don't have genes 
that are too related to each other.

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
Score, 5'-overhang, 3'-overhang, and Frame. Strand is +/-. Score is any 
floating point value. 5'- and 3'-overhang are the number of bp of an incomplete 
codon at each end of an exon. Frame is the reading frame (0..2 and *not* 1..3). 
Here's an example of the long format:

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


TUTORIAL
--------

In this tutorial, we will create SNAP HMM files for 3 different genomes. In the 
`DATA` directory, you will find fasta and gff3 files corresponding to 1 percent 
of the A. thaliana, C. elegans, and D. melanogaster genomes. Let's start by 
creating a directory for training A. thaliana in the main SNAP directory. We'll 
run `gff3_to_zff.pl` to convert the annotation to ZFF.

```
mkdir train_at
cd train_at
../gff3_to_zff.pl ../DATA/at.fa.gz ../DATA/at.gff3.gz > at.zff
```

The next step is to check for errors in the annotation. The training procedure
assumes that genes are _canonical_ in various respects.

+ Coding sequences start with ATG and end in a stop codon
+ Splice sites follow the GT..AG rule
+ Genes are at least 150 bp
+ Exons are at least 6 bp
+ CDS are at least 150 bp
+ Introns are at least 30 bp

Running `fathom -validate` will tell you which genes look ok and which genes 
look suspicious. Let's try one.

```
../fathom -validate at.zff ../DATA/at.fa.gz > at.validate
```

This will produce a bunch of output to STDERR. You will see several WARNING 
lines saying the DNA and annotation don't have the same definition lines. 
That's okay. The ZFF contains only the sequence id from the FASTA file and not 
the whole definition line present in the original FASTA. The last line gives 
some overall stats.

```
463 genes, 463 OK, 40 warnings, 0 errors
```

If you examine the `at.validate` file, you will see warnings for some short 
genes, short exons, non-canonical introns, etc. We won't be using these to 
train SNAP. To split genes into various categories use the `-categorize` 
function of `fathom` and give it a value for how much intergenic sequence you 
want on each side of a gene. For example `fathom -categorize 1000` attempts to 
put 1000 bp of genomic sequence on each side of a gene. However, if two genes 
are close to each other, say only 400 bp apart, they split the intergenic 
sequence and each get 200 bp of intergenic.

```
../fathom -categorize 100 at.zff ../DATA/at.fa.gz
```

This produces several new files:

+ `alt.*` contains genes with alternative splicing
+ `err.*` contains genes with errors
+ `olp.*` contains overlapping genes
+ `uni.*` contains unique genes
+ `wrn.*` contains genes with warnings

`fathom` doesn't want to create training files with alternative splicing. It 
could create a case of overtraining for those specific genes. If you have a lot 
of alternative splicing, you may want to remove all of the isoforms except for 
the main one. `fathom` also doesn't know what to do with overlapping genes 
because it requires genes to have intergenic sequence on either side. These 
tend to be rare. Genes with unusual features or outright errors are separated 
also.


The next step is to export all of the `uni` genes into their plus-stranded 
versions.

```
../fathom -export 100 -plus uni.*
```

This creates 4 new files:

+ `export.aa`  contains the amino acid sequences
+ `export.ann` contains the annotation in ZFF format
+ `export.dna` contains the sequences in FASTA format
+ `export.tx` contains the sequences of the coding sequences

Ideally, all of the proteins in `export.aa` start with `M` and end with `*`. 
Similarly, the `export.tx` files should start with `ATG` and end in a stop 
codon. All of the genes should validate without any reported warnings or 
errors.

```
../fathom -validate export.ann export.dna
```

The next step is to run `forge`, which will create a large number of model 
files.

```
../forge export.ann export.dna
```

Finally, run `hmm-assembler.pl` to glue the various models together to form an 
hmm parameter file. There are several options, but we'll just use the defaults.

```
./hmm-assembler.pl A.thaliana . > at.hmm
```

To verify this works, you can try it on the various fasta files we've used.

```
../snap at.hmm export.dna
../snap at.hmm uni.dna
../snap at.hmm ../DATA/at.fa.gz
```

---

Next, we are going to try a slightly different training procedure that might be 
better if you have a lot of genes that are getting stuffed into the `wrn.*` 
files. Let's rewind a bit.

```
../fathom -categorize 100 at.zff ../DATA/at.fa.gz
```

Let's see how many genes are in each category.

```
grep -c ">" *.ann
alt.ann:77
err.ann:0
olp.ann:0
uni.ann:236
wrn.ann:33
```

The standard procedure will have us training from the 236 genes in the `uni` 
category. To recover the genes in the `wrn` category, we'll just glue them to 
the `uni` and then export that for the training.

```
cat uni.ann wrn.ann > glue.ann
cat uni.dna wrn.dna > glue.dna
../fathom -export 100 -plus glue.ann glue.dna
../forge export.ann export.dna
../hmm-assembler.pl whatever .
```

---

What about all of the genes in the `alt` category? Those genes are reported to 
have multiple isoforms. Training from a gene with 10 isoforms would count that 
gene 10 times, so these are generally skipped. However, as more and more 
isoforms are found, this will become problematic. You will need some way to 
figure out which isoform is _canonical_ and delete the rest. I don't have an 
automated way to do that as each community has their own standards.

---

Try the training procedures for the C. elegans and D. melanogaster genomes 
next. Note that these training sets represent just 1% of each chromosome and 
are just for demonstration purposes.
