/*****************************************************************************\
memloop.c

test program to find memory leaks

Copyright (C) Ian Korf 2002-2013.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The MIT License (MIT) - opensource.org/licenses/MIT

\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zoe.h"

void loopVec (void);
void loopIVec (void);
void loopFVec (void);
void loopTVec (void);
void loopFeature (void);
void loopFeatureVec (void);
void loopHash (void);
void loopXtree (void);
void loopDNA (void);
void loopProtein (void);
void loopFastaFile (void);
void loopFeatureTable (const char *);
void loopCDS (void);
void loopDistribution (void);
void loopDuration (void);
void loopTransition (void);
void loopState (void);
void loopModel (const char *);
void loopScanner (const char *);
void loopCounter (const char *);

int Iterations = 1000;

int main (int argc, char *argv[]) {
	
	/* process commandline */
	if (argc == 1) {
		zoeO("usage: %s <test> [iterations - default 1000]\n", argv[0]);
		zoeM(stdout, 10,
			"tests:",
			"  -all",
			"  -Vec, -IVec, -FVec, -TVec, -Hash -Xtree",
			"  -Feature, -FeatureVec",
			"  -FastaFile, -DNA, -Protein, -CDS",
			"  -Distribution, -Duration, -State, -Transition",
			"The following tests are not part of -all and require files",
			"  -FeatureTable <file>",
			"  -Model, -Scanner, -Counter",
			"  -HMM <file>");
		exit(1);
	}
	
	zoeSetProgramName(argv[0]);
	zoeSetOption("-all", 0);
	
	/* zoeTools tests */
	zoeSetOption("-Vec", 0);
	zoeSetOption("-IVec", 0);
	zoeSetOption("-FVec", 0);
	zoeSetOption("-TVec", 0);
	zoeSetOption("-Hash", 0);
	zoeSetOption("-Xtree", 0);
	
	/* zoeFeature tests */
	zoeSetOption("-Feature", 0);
	zoeSetOption("-FeatureVec", 0);
	
	/* sequence tests */
	zoeSetOption("-DNA", 0);
	zoeSetOption("-Protein", 0);
	zoeSetOption("-FastaFile", 0);
	
	/* HMM component tests */
	zoeSetOption("-Distribution", 0);
	zoeSetOption("-Duration", 0);
	zoeSetOption("-State", 0);
	zoeSetOption("-Transition", 0);
	
	/* CDS */
	zoeSetOption("-CDS", 0);

	/* file-based tests */
	zoeSetOption("-FeatureTable", 1);
	zoeSetOption("-Model", 1);
	zoeSetOption("-Scanner", 1);
	zoeSetOption("-Counter", 1);
	
	/* incomplete tests */
	zoeSetOption("-HMM", 1);
	
	/* command line */
	zoeParseOptions(&argc, argv);
	if (argc == 2) Iterations = atoi(argv[1]);
	if (zoeOption("-Vec")          || zoeOption("-all")) loopVec();
	if (zoeOption("-IVec")         || zoeOption("-all")) loopIVec();
	if (zoeOption("-FVec")         || zoeOption("-all")) loopFVec();
	if (zoeOption("-TVec")         || zoeOption("-all")) loopTVec();
	if (zoeOption("-Hash")         || zoeOption("-all")) loopHash();
	if (zoeOption("-Xtree")        || zoeOption("-all")) loopXtree();
	if (zoeOption("-Feature")      || zoeOption("-all")) loopFeature();
	if (zoeOption("-FeatureVec")   || zoeOption("-all")) loopFeatureVec();
	if (zoeOption("-DNA")          || zoeOption("-all")) loopDNA();
	if (zoeOption("-Protein")      || zoeOption("-all")) loopProtein();
	if (zoeOption("-FastaFile")    || zoeOption("-all")) loopFastaFile();
	if (zoeOption("-Distribution") || zoeOption("-all")) loopDistribution();
	if (zoeOption("-Duration")     || zoeOption("-all")) loopDuration();
	if (zoeOption("-State")        || zoeOption("-all")) loopState();
	if (zoeOption("-Transition")   || zoeOption("-all")) loopTransition();
	if (zoeOption("-CDS")          || zoeOption("-all")) loopCDS();
	if (zoeOption("-FeatureTable"))   loopFeatureTable(zoeOption("-FeatureTable"));
	if (zoeOption("-Model"))          loopModel(zoeOption("-Model"));
	if (zoeOption("-Scanner"))        loopScanner(zoeOption("-Scanner"));
	if (zoeOption("-Counter"))        loopCounter(zoeOption("-Counter"));
	
	return 0;
}

void loopVec (void) {
	int    i, j;
	zoeVec vec = NULL;
	
	zoeE("Vec\n");
	for (j = 0; j < Iterations; j++) {
		vec = zoeNewVec();
		for (i = 0; i < 30000; i++) zoePushVec(vec, NULL);
		zoeDeleteVec(vec);
	}
}

void loopIVec (void) {
	int     i, j;
	zoeIVec vec = NULL;
	
	zoeE("IVec\n");
	for (j = 0; j < Iterations; j++) {
		vec = zoeNewIVec();
		for (i = 0; i < 30000; i++) zoePushIVec(vec, i);
		zoeDeleteIVec(vec);
	}
}

void loopFVec (void) {
	int     i, j;
	zoeFVec vec = NULL;
	
	zoeE("FVec\n");
	for (j = 0; j < Iterations; j++) {
		vec = zoeNewFVec();
		for (i = 0; i < 30000; i++) zoePushFVec(vec, (float)i);
		zoeDeleteFVec(vec);
	}
}

void loopTVec (void) {
	int     i, j;
	zoeTVec vec = NULL;
	
	zoeE("TVec\n");
	for (j = 0; j < Iterations; j++) {
		vec = zoeNewTVec();
		for (i = 0; i < 2500; i++) zoePushTVec(vec, "foobar");
		zoeDeleteTVec(vec);
	}
}

void loopFeature (void) {
	int        i, j;
	zoeFeature f = NULL;
	
	zoeE("Feature\n");
	for (j = 0; j < Iterations; j++) {
		for (i = 0; i < 1500; i++) {
			f = zoeNewFeature(Exon, 100, 200, '+', 100, 0, 0, 0, "foobar"/*, NULL*/);
			zoeDeleteFeature(f);
		}
	}
}

void loopFeatureVec (void) {
	int           i, j;
	zoeFeatureVec vec = NULL;
	zoeFeature    f   = zoeNewFeature(Exon, 100, 200, '+', 100, 0, 0, 0, "foobar");
	
	zoeE("FeatureVec\n");
	for (j = 0; j < Iterations; j++) {
		vec = zoeNewFeatureVec();
		for (i = 0; i < 1500; i++) zoePushFeatureVec(vec, f);
		zoeDeleteFeatureVec(vec);
	}
}

void loopDNA (void) {
	int    i, j;
	zoeDNA dna = NULL;
	
	zoeE("DNA\n");
	for (j = 0; j < Iterations; j++) {
		for (i = 0; i < 625; i++) {
			dna = zoeNewDNA("foo", "ATAGCTAAT");
			zoeDeleteDNA(dna);
		}
	}	
}

void loopProtein (void) {
	int        i, j;
	zoeProtein pro = NULL;
	
	zoeE("Protein\n");
	for (j = 0; j < Iterations; j++) {
		for (i = 0; i < 1250; i++) {
			pro = zoeNewProtein("foo", "IAN");
			zoeDeleteProtein(pro);
		}
	}	
}

void loopFastaFile (void) {
	int           i, j;
	zoeFastaFile fasta = NULL;
	
	zoeE("FastaFile\n");
	for (j = 0; j < Iterations; j++) {
		for (i = 0; i < 1000; i++) {
			fasta = zoeNewFastaFile("foo", "any_text_will_do");
			zoeDeleteFastaFile(fasta);
		}
	}	
}

void loopHash (void) {
	int     i, j;
	zoeHash hash;
	char    * string = zoeMalloc(3);
	
	zoeE("Hash\n");
	string[0] = 'A';
	for (j = 0; j < Iterations; j++) {
		hash = zoeNewHash();
		for (i = 0; i < 3000; i++) {
			string[1] = i;
			string[2] = '\0';
			zoeSetHash(hash, string, string);
		}
		zoeDeleteHash(hash);
	}	
}

void loopXtree (void) {
	int      i, j, k;
	float    f;
	char     seq[8];
	zoeXtree xt;
	
	zoeE("Xtree\n");
	seq[7] = '\0';
	for (k = 0; k < Iterations; k ++) {
		xt = zoeNewXtree();
		for (j = 0; j < 200; j++) {
			for (i = 0; i < 7; i++) {
				f = (float)rand() / (float)RAND_MAX;
				if      (f < 0.25) seq[i] = 'A';
				else if (f < 0.50) seq[i] = 'C';
				else if (f < 0.75) seq[i] = 'G';
				else               seq[i] = 'T';
			}
			zoeSetXtree(xt, seq, (void*)1);
		}
		zoeDeleteXtree(xt);
	}
}

void loopFeatureTable (const char * file) {
	int             i, j;
	zoeFeatureTable ft;
	
	zoeE("FeatureTable\n");
	for (j = 0; j < Iterations; j++) {
		for (i = 0; i < 25; i++) {
			ft = zoeGetFeatureTable(file);
			zoeDeleteFeatureTable(ft);
		}
	}
}

void loopCDS (void) {
	int           i, j;
	zoeCDS        cds   = NULL;
	zoeFeatureVec exons = NULL;
	zoeDNA        dna   = NULL;
	zoeFeature    exon  = NULL;
		
	zoeE("CDS\n");
	for (i = 0; i < 100; i++) {
		dna = zoeNewDNA(">foo", "nnnnnATGATAGCGgtnnnnnnnnnnagATAGCGAATTAAnnnnnnnnnn");
		exons = zoeNewFeatureVec();
		exon  = zoeNewFeature(Einit, 5, 13, '+', 100, 0, 0, 2, "gene"/*, NULL*/);
		zoePushFeatureVec(exons, exon);
		zoeDeleteFeature(exon);
		exon  = zoeNewFeature(Eterm, 28, 39, '+', 100, 0, 0, 1, "gene"/*, NULL*/);
		zoePushFeatureVec(exons, exon);
		zoeDeleteFeature(exon);
		
		for (j = 0; j < Iterations; j++) {
			cds = NULL;
			cds = zoeNewCDS("foo", dna, exons);
			zoeDeleteCDS(cds);
		}
		zoeDeleteDNA(dna);
		zoeDeleteFeatureVec(exons);
	}
}

void loopDistribution (void) {
	int             i, j;
	zoeDistribution d;
	float           f[5] = {0, 1, 2, 3, 4};
	
	zoeE("Distribution\n");
	for (i = 0; i < 1750; i++) {
		for (j = 0; j < Iterations; j++) {
			d = zoeNewDistribution(DEFINED, 0, 5, 1, f);
			zoeDeleteDistribution(d);
		}
	}
}

void loopDuration (void) {
	int             i, j;
	float           f[5] = {0, 1, 2, 3, 4};
	zoeDistribution D[2];
	zoeDuration     d = NULL;
	
	zoeE("Duration\n");
	
	for (i = 0; i < 600; i++) {
		for (j = 0; j < Iterations; j++) {
			D[0] = zoeNewDistribution(DEFINED,  0,  5, 5, f);
			D[1] = zoeNewDistribution(CONSTANT, 6, -1, 1, f);
			d    = zoeNewDuration(Exon, 2, D);
			zoeDeleteDuration(d);
		}
	}
}

void loopState (void) {
	int      i, j;
	zoeState si, se, ss;
	
	zoeE("State\n");
	for (i = 0; i < 1000; i++) {
		for (j = 0; j < Iterations; j++) {
			si = zoeNewState(INTERNAL, Intron, 0.5, 0.5, 0, 0, 1);
			se = zoeNewState(EXTERNAL, Exon,   0.0, 0.0, 0, 0, 1);
			ss = zoeNewState(SHUTTLE,  Repeat, 0.0, 0.0, 0, 0, 1);
			zoeDeleteState(si);
			zoeDeleteState(se);
			zoeDeleteState(ss);
		}
	}
}

void loopTransition (void) {
	int           i, j;
	zoeTransition t;
	
	zoeE("Transition\n");
	for (i = 0; i < 750; i++) {
		for (j = 0; j < Iterations; j++) {
			t = zoeNewTransition("Exon", "Intron", 0.5);
			zoeDeleteTransition(t);	
		}
	}
}

void loopModel (const char * file) {
	int      i, j;
	zoeModel m = NULL;
	
	zoeE("Model\n");
	for (i = 0; i < 35; i++) {
		for (j = 0; j < Iterations; j++) {
			m = zoeGetModel(file);
			zoeDeleteModel(m);
		}
	}
}

void loopScanner (const char * file) {
	int        i, j;
	zoeScanner s = NULL;
	zoeModel   m = zoeGetModel(file);
	zoeDNA     d = zoeNewDNA("foo", "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
	zoeDNA     a = zoeAntiDNA("bar", d);
	
	zoeE("Scanner\n");
	for (i = 0; i < 350; i++) {
		for (j = 0; j < Iterations; j++) {
			s = zoeNewScanner(d, a, m);
			zoeSetScannerScore(s, 10, 10);
			zoeDeleteScanner(s);
		}
	}
	zoeDeleteModel(m);
	zoeDeleteDNA(d);
	zoeDeleteDNA(a);
}

void loopCounter (const char * file) {
	int        i, j;
	zoeCounter c = NULL;
	zoeModel   m = zoeGetModel(file);
	zoeDNA     d = zoeNewDNA("foo", "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
	zoeDNA     a = zoeAntiDNA("bar", d);
	
	zoeE("Counter\n");
	for (i = 0; i < 350; i++) {
		for (j = 0; j < Iterations; j++) {
			c = zoeNewCounter(d, a, m);
			zoeDeleteCounter(c);
		}
	}
	zoeDeleteModel(m);
	zoeDeleteDNA(d);
	zoeDeleteDNA(a);
}

