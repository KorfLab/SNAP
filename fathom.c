/*****************************************************************************\
 fathom.c

Annotation investigation utility

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
#include <time.h>
#include "zoe.h"

/* globals */
char * DNA_file   = NULL;
char * ANN_file   = NULL;
char * PRE_file   = NULL;
zoeFile DNA_stream;
zoeFile ANN_stream;
zoeFile PRE_stream;
zoeDNA          DNA   = NULL; /* the current DNA */
zoeDNA          ANTI  = NULL; /* reverse complement of DNA */
zoeFeatureTable ANN   = NULL; /* the current annotation */
zoeFeatureTable PRE   = NULL; /* the current predictions */
int PREDICTIONS   = 0; /* flag to determine if predictions are present */
int ANTI_REQUIRED = 0; /* flag to determine if ANTI DNA is required */

struct score_s {
	score_t   plus_score;
	score_t   minus_score;
	strand_t  strand;      /* +, -, =, 0 : 0 means pseudo-site on both strands */
};
typedef struct score_s score_s;

/* prototypes - showing relationships */
void help (void);
void anti (void);
void fuse (void);
void validate (void);
void reformat (void);
void genestats (void);
void split (void);
void exon_intron (void);
void extract ();
	zoeFeatureVec extractAcceptor (void);
	zoeFeatureVec extractDonor (void);
	zoeFeatureVec extractStart (void);
	zoeFeatureVec extractStop (void);
	zoeFeatureVec extractPolyA (void);
	zoeFeatureVec extractInter (void);
	zoeFeatureVec extractIntron (void);
	zoeFeatureVec extractUTR5 (void);
	zoeFeatureVec extractUTR3 (void);
	zoeFeatureVec extractExon (void);
	zoeFeatureVec extractFeature (zoeLabel);
void export (void);
	zoeHash kill_list (const zoeVec);
void categorize (void);
void export_feature (void);
void score_genes (void);
void filter_genes (void);
void compare_genes (void);
void ace_format (void);
void openData (void);
	void   closeData (void);
	int    getData (void);

static float gc_content (const zoeDNA dna) {
	int i;
	int count[5];
	
	for (i = 0; i < 5; i++) count[i] = 0;
	for (i = 0; i < dna->length; i++) count[(int)dna->s5[i]]++;
	return (float)(count[1]+count[2]) / (float)(count[0]+count[1]+count[2]+count[3]);
	
}

static char usage[]  = "\n\
FATHOM - sequence and annotation tool (version 2006-07-28)\n\n\
usage: fathom <ann> <dna> <commands>\n\
commands:\n\
  -help           report useful information\n\
  -validate [-quiet]\n\
  -gene-stats [-errors-ok -nucleotide -dinucleotide]\n\
  -categorize <int>\n\
  -export <int> [-plus -errors-ok]\n\
  -extract <feature> -length <int> -offset <int>\n\
  -exon-intron\n\
  -split <-number <int> | -training <float> | -GC <float> | -repeat <float>>\n\
  -ace-format <-gene-method <string> [-dna -extra <string>]>\n\
  -compare-genes <predictions> [-details]\n\
  -score-genes <hmm> [-errors-ok]\n\
  -filter-genes <hmm> -min-score <float> -min-length <int>\n\
";

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char **argv) {
	
	/* command line */
	
	zoeSetProgramName(argv[0]);
	
	/* dataset */
	zoeSetOption("-help",       0);
	zoeSetOption("-randseed",   1);
	zoeSetOption("-lcfilter",   0);
	zoeSetOption("-lcsmooth",   0);
	zoeSetOption("-validate",   0);
		zoeSetOption("-quiet",  0);
	zoeSetOption("-gene-stats", 0);
		zoeSetOption("-nucleotide",   0);
		zoeSetOption("-dinucleotide", 0);
	zoeSetOption("-reformat",   0);
	zoeSetOption("-ace-format",      0);
		zoeSetOption("-gene-method", 1);
		zoeSetOption("-extra",       1);
		zoeSetOption("-dna",         0);
	zoeSetOption("-extract",    1);
		zoeSetOption("-length", 1);
		zoeSetOption("-offset", 1);
	zoeSetOption("-export",             1);
		zoeSetOption("-plus",           0);
		zoeSetOption("-errors-ok",      0);
	zoeSetOption("-exon-intron", 0);
	zoeSetOption("-split",        0);
		zoeSetOption("-number",   1);
		zoeSetOption("-training", 1);
		zoeSetOption("-GC",       1);
		zoeSetOption("-repeat",   1);
	zoeSetOption("-anti", 0);
	zoeSetOption("-fuse", 1);
	zoeSetOption("-categorize",     1);
		zoeSetOption("-min-intron", 1);
		zoeSetOption("-max-intron", 1);
		zoeSetOption("-min-exon",   1);
		zoeSetOption("-max-exon",   1);
		zoeSetOption("-min-gene",   1);
		zoeSetOption("-max-gene",   1);
		zoeSetOption("-min-cds",    1);
		zoeSetOption("-max-cds",    1);
	
	/* scoring, comparison, filtering */
	zoeSetOption("-compare-genes",  1);
		zoeSetOption("-details",    0);
	zoeSetOption("-score-genes",    1);
	zoeSetOption("-filter-genes",   1);
		zoeSetOption("-min-score",  1);
		zoeSetOption("-min-length", 1);
	
	/* begin */
	zoeParseOptions(&argc, argv);
	if (zoeOption("-help")) help();
	
	if (argc != 3) {
		zoeE("%s", usage);
		exit(1);
	}
	
	ANN_file = argv[1];
	DNA_file = argv[2];
	
	if (zoeOption("-compare-genes")) {
		PRE_file = zoeOption("-compare-genes");
		PREDICTIONS = 1;
	}
	
	if (zoeOption("-score-genes")
		|| zoeOption("-anti")
		|| zoeOption("-fuse")
		|| zoeOption("-filter-genes")
	) ANTI_REQUIRED = 1;
	
	if (zoeOption("-filter-genes")) {
		if (!zoeOption("-min-score"))  zoeExit("-min-score required");
		if (!zoeOption("-min-length")) zoeExit("-min-length required");
	}
	
	if (zoeOption("-randseed")) srand(atoi(zoeOption("-randseed")));
	else                        srand(time(NULL));
		
	if (zoeOption("-min-intron")) zoeSetMinIntron(atoi(zoeOption("-min-intron")));	
	if (zoeOption("-max-intron")) zoeSetMaxIntron(atoi(zoeOption("-max-intron")));
	if (zoeOption("-min-exon"))   zoeSetMinExon(atoi(zoeOption("-min-exon")));	
	if (zoeOption("-max-exon"))   zoeSetMaxExon(atoi(zoeOption("-max-exon")));
	if (zoeOption("-min-gene"))   zoeSetMinGene(atoi(zoeOption("-min-gene")));	
	if (zoeOption("-max-gene"))   zoeSetMaxGene(atoi(zoeOption("-max-gene")));
	if (zoeOption("-min-cds"))    zoeSetMinCDS(atoi(zoeOption("-min-cds")));	
	if (zoeOption("-max-cds"))    zoeSetMaxCDS(atoi(zoeOption("-max-cds")));
	
	/* main control */
	if      (zoeOption("-validate"))      validate();
	else if (zoeOption("-gene-stats"))    genestats();
	else if (zoeOption("-reformat"))      reformat();
	else if (zoeOption("-categorize"))    categorize();
	else if (zoeOption("-export"))        export();
	else if (zoeOption("-exon-intron"))   exon_intron();
	else if (zoeOption("-anti"))          anti();
	else if (zoeOption("-fuse"))          fuse();
	else if (zoeOption("-score-genes"))   score_genes();
	else if (zoeOption("-filter-genes"))  filter_genes();
	else if (zoeOption("-compare-genes")) compare_genes();
	else if (zoeOption("-split"))         split();	
	else if (zoeOption("-ace-format"))    ace_format();
	else if (zoeOption("-extract"))       extract();
		
	return 0;
}

void anti (void) {
	FILE *DNA_out, *ANN_out;
		
	DNA_out = fopen("anti.dna", "w");
	ANN_out = fopen("anti.ann", "w");
		
	openData();
	while (getData()) {
		zoeAntiFeatureTable(ANN, DNA->length);
		zoeWriteTriteFeatureTable(ANN_out, ANN);
		zoeWriteDNA(DNA_out, ANTI);
	}
	closeData();
	
	fclose(DNA_out);
	fclose(ANN_out);
}

void fuse (void) {
	int i, j, group_size, group_num, offset;
	float f;
	FILE *DNA_out, *ANN_out;
	
	DNA_out = fopen("fuse.dna", "w");
	ANN_out = fopen("fuse.ann", "w");
	
	group_size = atoi(zoeOption("-fuse"));
	if (group_size < 1) zoeExit("-fuse must be >= 1");
	
	group_num = 0;
	openData();
	while (getData()) {
		offset = 0;
		group_num++;
		zoeS(DNA_out, ">fuse-%d\n", group_num);
		zoeS(ANN_out, ">fuse-%d\n", group_num);
		
			f = (float)rand()/(float)RAND_MAX;
			if (f >= 0.5) {
				zoeS(DNA_out, "%s\n", DNA->seq);
			} else {
				zoeS(DNA_out, "%s\n", ANTI->seq);
				zoeAntiFeatureTable(ANN, DNA->length);
			}
			for (j = 0; j < ANN->vec->size; j++) {
				ANN->vec->elem[j]->start += offset;
				ANN->vec->elem[j]->end   += offset;
				zoeWriteTriteFeature(ANN_out, ANN->vec->elem[j]);
			}
			offset += DNA->length;
		
		for (i = 1; i < group_size; i++) {
			if (getData() == 0) break;
			
			f = (float)rand()/(float)RAND_MAX;
			if (f >= 0.5) {
				zoeS(DNA_out, "%s\n", DNA->seq);
			} else {
				zoeS(DNA_out, "%s\n", ANTI->seq);
				zoeAntiFeatureTable(ANN, DNA->length);
			}
			for (j = 0; j < ANN->vec->size; j++) {
				ANN->vec->elem[j]->start += offset;
				ANN->vec->elem[j]->end   += offset;
				zoeWriteTriteFeature(ANN_out, ANN->vec->elem[j]);
			}
			
			offset += DNA->length;
		}
	}
	closeData();
	
	fclose(DNA_out);
	fclose(ANN_out);

}

/*****************************************************************************\
 validate
\*****************************************************************************/

void validate (void) {
	int       i, tested = 0, OK = 0, warnings = 0, errors = 0;
	zoeCDS    cds;
	zoeVec    genes;
	int       quiet = zoeOption("-quiet") ? 1 : 0;
	
	openData();
	while (getData()) {
		if (!quiet && strcmp(DNA->def, ANN->def) != 0) {
			zoeWarn("DNA and annotation do not have the same def line!\n");
			zoeWarn("%s\n%d\n", DNA->def, ANN->def);
		}
						
		genes = zoeGetGenes(ANN, DNA);
		for (i = 0; i < genes->size; i++) {
			tested++;
			cds = genes->elem[i];
			if (cds->errors->size)   errors++;
			if (cds->warnings->size) warnings++;
			if (cds->OK)             OK++;
			
			zoeO("%s: ", DNA->def);
			if (cds->errors->size || cds->warnings->size) {
				zoeReportCDS(stdout, cds);
				zoeWriteCDS(stdout, cds);
			} else {
				zoeO("%s OK\n", cds->name);
			}
			zoeDeleteCDS(cds);
			
		}
		zoeDeleteVec(genes);
	}
	closeData();
	zoeE("%d genes, %d OK, %d warnings, %d errors\n", tested, OK, warnings, errors);
}


/*****************************************************************************\
 reformat
\*****************************************************************************/

void reformat (void) {
	int               i;
	zoeFeatureTable   ft;
	zoeCDS            cds;
	zoeVec            genes;

	
	openData();
	while (getData()) {
		genes = zoeGetGenes(ANN, DNA);
		for (i = 0; i < genes->size; i++) {
			cds = genes->elem[i];
			if (cds->OK) {
				if (zoeOption("-canonical") && cds->warnings->size) {
					zoeE("%s skipped due to warnings");
				} else {
					ft = zoeNewFeatureTable(cds->name, cds->exons);
					zoeWriteFeatureTable(stdout, ft);
					zoeDeleteFeatureTable(ft);
				}
			} else {
				zoeE("%s skipped due to errors\n", cds->name);
			}
			zoeDeleteCDS(cds);
		}
		zoeDeleteVec(genes);
	}
	closeData();
}

/*****************************************************************************\
 gene stats
\*****************************************************************************/

float entropy (zoeDNA dna, zoeFeature exon) {
	int        i;
	float      f, h, H = 0;
	int        symbol[256];
	zoeProtein pro = zoeTranslateFeature(">foo", dna, exon);
	
	
	for (i = 0; i < 256; i++) symbol[i] = 0;
	for (i = 0; i < pro->length; i++) symbol[(int)pro->seq[i]]++;
	for (i = 0; i < 256; i++) {
		if (symbol[i]) {			
			f = (float)symbol[i] / (float)pro->length;
			h = f * log(f);
			H -= h;
		}
	}
	
	zoeDeleteProtein(pro);
	return H;
}

void genestats (void) {
	int     i, j, length;
	int     gene_cnt = 0, max_gene = 0, min_gene = INT_MAX, plus = 0, minus = 0;
	int     Esngl = 0, multi = 0;
	int     exons   = 0, exon_len   = 0, max_exon   = 0, min_exon   = INT_MAX;
	int     introns = 0, intron_len = 0, max_intron = 0, min_intron = INT_MAX;
	int     dinuc[25];
	int     nuc[5];
	zoeCDS  gene;
	zoeVec  genes;
	zoeFeature exon, intron;
	int     seq_count = 0;
	float   h, gc, min_gc = 1, max_gc = 0, gc_total = 0;
	int     H[50];
	int     index, total;
	char *  nt[5] = {"A", "C", "G", "T", "N"};
	char *  dd[25] = {
		"AA", "AC", "AG", "AT", "AN",
		"CA", "CC", "CG", "CT", "CN",
		"GA", "GC", "GG", "GT", "GN",
		"TA", "TC", "TG", "TT", "TN",
		"NA", "NC", "NG", "NT", "NN"};
	
	for (i = 0; i < 5; i++) nuc[i] = 0;
	for (i = 0; i < 25; i++) dinuc[i] = 0;
	for (i = 0; i < 50; i++) H[i] = 0;

	openData();
	while (getData()) {
		if (zoeOption("-nucleotide")) {
			for (i = 0; i < DNA->length; i++) {
				nuc[(int)DNA->s5[i]]++;
			}
		}
		
		if (zoeOption("-dinucleotide")) {
			for (i = 1; i < DNA->length; i++) {
				index = DNA->s5[i] + (5 * DNA->s5[i-1]);
				dinuc[index]++;
			}
		}
	
		gc = gc_content(DNA);
		if (gc > max_gc) max_gc = gc;
		if (gc < min_gc) min_gc = gc;	
		gc_total += gc;
		seq_count++;
	
		genes = zoeGetGenes(ANN, DNA);
		for (i = 0; i < genes->size; i++) {
			gene = genes->elem[i];
			
			if (!gene->OK && !zoeOption("-errors-ok")) {
				zoeE("%s skipped due to errors\n", gene->name);
				continue;
			} else if (zoeOption("-canonical") && gene->warnings->size) {
				zoeE("%s skipped due to warnings\n", gene->name);
				continue;
			}
						
			length = gene->end - gene->start + 1;
			if (length < min_gene) min_gene = length;
			if (length > max_gene) max_gene = length;
						
			gene_cnt++;
			
			if (gene->strand == '+') plus++;
			else                     minus++;
			
			if (gene->exons->size == 1) Esngl++;
			else                        multi++;
			
			exons += gene->exons->size;
			introns += gene->introns->size;
			
			for (j = 0; j < gene->exons->size; j++) {
				exon = gene->exons->elem[j];
				length = exon->end - exon->start + 1;
				exon_len += length;
				if (length < min_exon) min_exon = length;
				if (length > max_exon) max_exon = length;
				if (zoeOption("-entropy")) {
					h = entropy(gene->dna, exon);
					H[ (int)(h * 10) ]++;
				}
			}
			
			for (j = 0; j < gene->introns->size; j++) {
				intron = gene->introns->elem[j];
				length = intron->end - intron->start +1;
				intron_len += length;
				if (length < min_intron) min_intron = length;
				if (length > max_intron) max_intron = length;
			}
		}	
		for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
		zoeDeleteVec(genes);
	}
	closeData();
	
	zoeO("%d sequences\n", seq_count);
	zoeO("%f avg GC fraction (min=%f max=%f)\n", gc_total/(float)seq_count, min_gc, max_gc);
	zoeO("%d genes (plus=%d minus=%d)\n", gene_cnt, plus, minus);
	zoeO("%d (%f) single-exon\n", Esngl, (float)Esngl/(float)gene_cnt);
	zoeO("%d (%f) multi-exon\n", multi, (float)multi/(float)gene_cnt);
	zoeO("%f mean exon (min=%d max=%d)\n", (float)exon_len / (float)exons, min_exon, max_exon);
	zoeO("%f mean intron (min=%d max=%d)\n", (float)intron_len / (float)introns, min_intron, max_intron);
	
	if (zoeOption("-entropy")) {
		zoeO("Exon entropy\n");
		for (i = 0; i < 50; i++) {
			if (H[i]) zoeO("%.1f %d\n", (float)i/10, H[i]);
		}
	}
	
	if (zoeOption("-nucleotide")) {
		total = 0;
		for (i = 0; i < 5; i++) total += nuc[i];
		for (i = 0; i < 5; i++) {
			zoeO("%s\t%f\n", nt[i], (float)nuc[i] / (float)total);
		}
	}
	if (zoeOption("-dinucleotide")) {	
		total = 0;
		for (i = 0; i < 25; i++) total += dinuc[i];
		for (i = 0; i < 25; i++) {
			zoeO("%s\t%f\n", dd[i], (float)dinuc[i] / (float)total);
		}
	}
}

/*****************************************************************************\
 exon-intron
\*****************************************************************************/

void exon_intron (void) {
	int    i, j;
	zoeCDS cds;
	zoeVec genes;
	zoeDNA dna;
	
	openData();
	while (getData()) {
		genes = zoeGetGenes(ANN, DNA);
		if (genes->size != 1) {
			zoeExit("fathom -exon-intron requires 1 gene per sequence");
		}

		cds = genes->elem[0];
		zoeO(">%s\n", cds->name);
		for (i = 0; i < cds->exons->size; i++) {
			dna = zoeFeatureDNA("exon", cds->dna, cds->exons->elem[i]);
			for (j = 0; j < dna->length; j++) {
				if (dna->seq[j] > 'Z') dna->seq[j] -= 32;
			}
			zoeO("%s\n", dna->seq);
			zoeDeleteDNA(dna);
			
			if (cds->introns->size && i < cds->introns->size) {
				dna = zoeFeatureDNA("intron", cds->dna, cds->introns->elem[i]);
				for (j = 0; j < dna->length; j++) {
					if (dna->seq[j] < 'a') dna->seq[j] += 32;
				}
				zoeO("%s\n", dna->seq);
				zoeDeleteDNA(dna);
			}
		}
		zoeO("\n");
		zoeDeleteCDS(cds);
		zoeDeleteVec(genes);
	}
	closeData();
	
}


/*****************************************************************************\
 split functions
\*****************************************************************************/

float gc_frac (void) {
	int i, c[5];
	
	for (i = 0; i < 5; i++) c[i] = 0;
	for (i = 0; i < DNA->length; i++) c[(int)DNA->s5[i]]++;
	
	return (float)(c[1] + c[2]) / (float)(c[0] + c[1] + c[2] + c[3]);
}

void split_by_number (void) {
	FILE *dna[20], *ann[20];
	int  split;
	int  i;
	char filename[256];
	int  set;
	
	split = atoi(zoeOption("-number"));
	if (split < 2) {
		zoeE("setting split -number to minimum value of 2\n");
		split = 2;
	} else if (split > 20) {
		zoeE("setting split -number to maximum value of 20\n");
		split = 20;
	}
	
	for (i = 0; i < split; i++) {
		sprintf(filename, "%d.dna", i);
		if ((dna[i] = fopen(filename, "w")) == NULL) zoeExit("error opening %s", filename);
		sprintf(filename, "%d.ann", i);
		if ((ann[i] = fopen(filename, "w")) == NULL) zoeExit("error opening %s", filename);
	}
	
	i = 0;
	openData();
	while (getData()) {
		set = i % split;
		zoeWriteDNA(dna[set], DNA);
		zoeWriteFeatureTable(ann[set], ANN);
		i++;
	}
	closeData();
	
	for (i = 0; i < split; i++) {
		fclose(dna[i]);
		fclose(ann[i]);
	}
	
}

void split (void) {
	float f, cutoff;
	FILE *train_dna, *train_ann, *test_dna, *test_ann;
	FILE *hi_dna, *hi_ann, *lo_dna, *lo_ann;
	
	train_dna = train_ann = test_dna = test_ann = hi_dna = hi_ann = lo_dna = lo_ann = NULL;
	
	if (zoeOption("-number")) {
		split_by_number();
		return;
	}
	
	if      (zoeOption("-training")) cutoff = atof(zoeOption("-training"));
	else if (zoeOption("-GC"))       cutoff = atof(zoeOption("-GC"));
	else {
		cutoff = 0; /* shush */
		zoeExit("no data splitting method given, try -GC or -training");
	}
	if (cutoff <= 0 || cutoff >= 1) zoeExit("-split (%f) is out of range", cutoff);

	/* open output files as "training" and "testing" */
	if (zoeOption("-training")) {
		if ((train_dna = fopen("train.dna", "w")) == NULL) zoeExit("can't open train.dna");
		if ((train_ann = fopen("train.ann", "w")) == NULL) zoeExit("can't open train.ann");
		if ((test_dna  = fopen("test.dna",  "w")) == NULL) zoeExit("can't open test.dna");
		if ((test_ann  = fopen("test.ann",  "w")) == NULL) zoeExit("can't open test.ann");
	} else {
		if ((hi_dna = fopen("hi-gc.dna", "w")) == NULL) zoeExit("can't open hi-gc.dna");
		if ((hi_ann = fopen("hi-gc.ann", "w")) == NULL) zoeExit("can't open hi-gc.ann");
		if ((lo_dna = fopen("lo-gc.dna", "w")) == NULL) zoeExit("can't open lo-gc.dna");
		if ((lo_ann = fopen("lo-gc.ann", "w")) == NULL) zoeExit("can't open lo-gc.ann");
	}
	
	/* repeat content not yet implemented */

	openData();
	while (getData()) {
		if      (zoeOption("-training")) f = (float)rand()/(float)RAND_MAX;
		else if (zoeOption("-GC"))       f = gc_frac();
		else                             f = 0; /* shush */
		
		if (f < cutoff) {
			if (zoeOption("-training")) {
				zoeWriteDNA(train_dna, DNA);
				zoeWriteTriteFeatureTable(train_ann, ANN);
			} else {
				zoeWriteDNA(lo_dna, DNA);
				zoeWriteTriteFeatureTable(lo_ann, ANN);
			}
		} else {
			if (zoeOption("-training")) {
				fflush(test_ann);
				zoeWriteTriteFeatureTable(test_ann, ANN);
				zoeWriteDNA(test_dna, DNA);
			} else {
				zoeWriteDNA(hi_dna, DNA);
				zoeWriteTriteFeatureTable(hi_ann, ANN);
			}
		}
	}
	closeData();
	
	if (zoeOption("-training")) {
		fclose(train_dna);
		fclose(train_ann);
		fclose(test_dna);
		fclose(test_ann);
	} else {
		fclose(hi_dna);
		fclose(hi_ann);
		fclose(lo_dna);
		fclose(lo_ann);
	}
}

/*****************************************************************************\
 extract
\*****************************************************************************/

void extract (void) {
	int               i, length, offset, point_feature, id;
	char            * name, token[64];
	zoeLabel          label;
	zoeFeatureVec     vec = NULL;
	zoeFeature        f;
	zoeDNA            seq;
	zoeVec            genes;
	zoeCDS            gene;
	
	name = zoeOption("-extract");
	label = zoeText2Label(zoeOption("-extract"));
	
	/* point features requires length and offset */
	switch (label) {
		case Acceptor: case Donor: case Start: case Stop:
			if (! (zoeOption("-length") && zoeOption("-offset")))
				zoeExit("extract requires -length and -offset for %s", name);
			length = atoi(zoeOption("-length"));
			offset = atoi(zoeOption("-offset"));
			point_feature = 1;
			break;
		case Coding: case Intron: case Inter:
			if (zoeOption("-length") || zoeOption("-offset"))
				zoeWarn("ignoring -length and -offset for %s", name);
			length = 0; offset = 0; point_feature = 0; /* shush */
			break;
		default:
			length = 0; offset = 0; point_feature = 0; /* shush */
			zoeExit("extract does not support %s yet", name);
	}

	
	/* main loop */
	id = 0; /* index number for features */
	openData();
	while (getData()) {
		switch (label) {
			case Acceptor: vec = extractAcceptor();     break;
			case Donor:    vec = extractDonor();        break;
			case Start:    vec = extractStart();        break;
			case Stop:     vec = extractStop();         break;
			case PolyA:    vec = extractFeature(PolyA); break;
			case Inter:    vec = extractInter();        break;
			case Intron:   vec = extractIntron();       break;
			case UTR5:     vec = extractUTR5();         break;
			case UTR3:     vec = extractUTR3();         break;
			case Coding:
				genes = zoeGetGenes(ANN, DNA);
				for (i = 0; i < genes->size; i++) {
					gene = genes->elem[i];
					zoeWriteDNA(stdout, gene->tx);
					zoeDeleteCDS(gene);
				}
				zoeDeleteVec(genes);
				continue; /* back to while */
			default:
				zoeExit("odd extract error"); 
		}
		
		for (i = 0; i < vec->size; i++) {
			f = vec->elem[i];
			if (f->start - offset < 0) continue;
			
			id++;
			sprintf(token, "%s-%d from sequence %s", name, id, DNA->def);
			
			if (point_feature) seq = zoeSubseqDNA(token, DNA, f->start - offset, length);
			else               seq = zoeFeatureDNA(token, DNA, f);
			zoeWriteDNA(stdout, seq);
			zoeDeleteDNA(seq);
		}
		zoeDeleteFeatureVec(vec);
	}
	
	closeData();
}

zoeFeatureVec extractInter (void) {
	int           i, start, end;
	zoeCDS        gene, prev_gene, this_gene;
	zoeFeature    inter_p, inter_m;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	/* right now, this function ignores UTRs */
		
	/* space before first gene */
	gene = genes->elem[0];
	start = 0;
	end   = gene->start -1;
	if (end >= start) {
		inter_p = zoeNewFeature(Inter, start, end, '+', MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, "before"/*, NULL*/);
		inter_m = zoeNewFeature(Inter, start, end, '-', MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, "before"/*, NULL*/);
		zoePushFeatureVec(vec, inter_p);
		zoePushFeatureVec(vec, inter_m);
		zoeDeleteFeature(inter_p);
		zoeDeleteFeature(inter_m);
	}
	
	/* between genes */
	for (i = 1; i < genes->size; i++) {
		prev_gene = genes->elem[i-1];
		this_gene = genes->elem[i];
		start     = prev_gene->end +1;
		end       = this_gene->start -1;
		if (end >= start) {
			inter_p = zoeNewFeature(Inter, start, end, '+', MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, "between"/*, NULL*/);
			inter_m = zoeNewFeature(Inter, start, end, '-', MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, "between"/*, NULL*/);
			zoePushFeatureVec(vec, inter_p);
			zoePushFeatureVec(vec, inter_m);
			zoeDeleteFeature(inter_p);
			zoeDeleteFeature(inter_m);
		}
	}
	
	/* space after genes */
	gene = genes->elem[genes->size -1];
	start = gene->end +1;
	end   = DNA->length -1;
	if (end >= start) {
		inter_p = zoeNewFeature(Inter, start, end, '+', MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, "after"/*, NULL*/);
		inter_m = zoeNewFeature(Inter, start, end, '-', MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, "after"/*, NULL*/);
		zoePushFeatureVec(vec, inter_p);
		zoePushFeatureVec(vec, inter_m);
		zoeDeleteFeature(inter_p);
		zoeDeleteFeature(inter_m);
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractUTR5 (void) {
	int           i, offset = 0, length = 0, start, end;
	zoeCDS        gene;
	zoeFeature    utr5;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	if (zoeOption("-offset")) offset = atoi(zoeOption("-offset"));
	else zoeExit("you must supply -offset with UTR5/Upstream");
	
	if (zoeOption("-length")) length = atoi(zoeOption("-length"));
	else zoeExit("you must supply -length with UTR5/Upstream");
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		if (!gene->end_found) continue;
		if (gene->strand == '+') {
			end   = gene->start -1 -offset;
			start = end - length;
			if (end < 0) continue; /* no UTR5 possible */
			if (start < 0) start = 0;
		} else {
			start = gene->end +1 + offset;
			end   = start + length;
			if (start > DNA->length -1) continue; /* no UTR5 possible */
			if (end > DNA->length -1) end = DNA->length -1;
		}
		utr5 = zoeNewFeature(UTR5, start, end, gene->strand, MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, NULL/*, NULL*/);
		zoePushFeatureVec(vec, utr5);
		zoeDeleteFeature(utr5);
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractUTR3 (void) {
	int           i, offset = 0, length = 0, start, end;
	zoeCDS        gene;
	zoeFeature    utr3;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	if (zoeOption("-offset")) offset = atoi(zoeOption("-offset"));
	else zoeExit("you must supply -offset with UTR3/Downstream");
	
	if (zoeOption("-length")) length = atoi(zoeOption("-length"));
	else zoeExit("you must supply -length with UTR3/Downstream");
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		if (!gene->end_found) continue;
		if (gene->strand == '+') {
			start = gene->end +1 + offset;
			end   = start + length;
			if (start > DNA->length -1) continue; /* no UTR3 possible */
			if (end > DNA->length -1) end = DNA->length -1;
		} else {
			end   = gene->start -1 -offset;
			start = end - length;
			if (end < 0) continue; /* no UTR3 possible */
			if (start < 0) start = 0;
		}
		utr3 = zoeNewFeature(UTR3, start, end, gene->strand, MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, NULL/*, NULL*/);
		zoePushFeatureVec(vec, utr3);
		zoeDeleteFeature(utr3);
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractAcceptor (void) {
	int           i, j, coor;
	zoeCDS        gene;
	zoeFeature    acc, exon;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->exons->size; j++) {
			exon = gene->exons->elem[j];
			if ( !(exon->label == Exon || exon->label == Eterm) ) continue;
			if (exon->strand == '+') coor = exon->start -1;
			else                     coor = exon->end +1;
			acc = zoeNewFeature(Acceptor, coor, coor, exon->strand, MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, gene->name/*, NULL*/);
			zoePushFeatureVec(vec, acc);
			zoeDeleteFeature(acc);
		}
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractDonor (void) {
	int           i, j, coor;
	zoeCDS        gene;
	zoeFeature    donor, exon;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->exons->size; j++) {
			exon = gene->exons->elem[j];
			if ( !(exon->label == Exon || exon->label == Einit) ) continue;
			if (exon->strand == '+') coor = exon->end +1;
			else                     coor = exon->start -1;
			donor = zoeNewFeature(Donor, coor, coor, exon->strand, MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, gene->name/*, NULL*/);
			zoePushFeatureVec(vec, donor);
			zoeDeleteFeature(donor);
		}
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractStart (void) {
	int           i, j, coor;
	zoeCDS        gene;
	zoeFeature    start, exon;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->exons->size; j++) {
			exon = gene->exons->elem[j];
			if ( !(exon->label == Einit || exon->label == Esngl) ) continue;
			if (exon->strand == '+') coor = exon->start;
			else                     coor = exon->end;
			start = zoeNewFeature(Start, coor, coor, exon->strand, MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, gene->name/*, NULL*/);
			zoePushFeatureVec(vec, start);
			zoeDeleteFeature(start);
		}
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractStop (void) {
	int           i, j, coor;
	zoeCDS        gene;
	zoeFeature    stop, exon;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->exons->size; j++) {
			exon = gene->exons->elem[j];
			if ( !(exon->label == Eterm || exon->label == Esngl) ) continue;
			if (exon->strand == '+') coor = exon->end -2;
			else                     coor = exon->start +2;
			stop = zoeNewFeature(Stop, coor, coor, exon->strand, MIN_SCORE, UNDEFINED_FRAME, UNDEFINED_FRAME, UNDEFINED_FRAME, gene->name/*, NULL*/);
			zoePushFeatureVec(vec, stop);
			zoeDeleteFeature(stop);
		}
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractPolyA (void) {
	int           i;
	zoeFeature    f;
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	for (i = 0; i < ANN->vec->size; i++) {
		f = ANN->vec->elem[i];
		if (f->label == PolyA) zoePushFeatureVec(vec, f);
	}
	
	return vec;
}

zoeFeatureVec extractIntron (void) {
	int           i, j;
	zoeCDS        gene;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	zoeFeatureVec frags;
	zoeFeature    intron, tmp;
	int           start, end;
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->introns->size; j++) {
			zoePushFeatureVec(vec, gene->introns->elem[j]);
		}
	}
	
	/* fragment Ns? */
	if (zoeOption("-fragmentN")) {
		frags = zoeNewFeatureVec();
		for (i = 0; i < vec->size; i++) {
			intron = vec->elem[i];
			start = intron->start;
			end   = start;
			while (end < intron->end) {
				if (DNA->s5[end] != 4) end++;
				else {
					tmp = zoeNewFeature(Intron, start, end, '+', 0, 0, 0, 0, NULL/*, NULL*/);
					zoePushFeatureVec(frags, tmp);
					zoeDeleteFeature(tmp);
					start = end;
					while (DNA->s5[start] == 4) start++;
					end = start;
				}
			}
			tmp = zoeNewFeature(Intron, start, end, '+', 0, 0, 0, 0, NULL/*, NULL*/);
			zoePushFeatureVec(frags, tmp);
			zoeDeleteFeature(tmp);
		}
		
		zoeDeleteFeatureVec(vec);
		vec = frags;
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractExon (void) {
	int           i, j;
	zoeCDS        gene;
	zoeVec        genes = zoeGetGenes(ANN, DNA);
	zoeFeatureVec vec   = zoeNewFeatureVec();
	
	/* collect features */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->exons->size; j++) {
			zoePushFeatureVec(vec, gene->exons->elem[j]);
		}
	}
	
	/* end */
	for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
	zoeDeleteVec(genes);
	return vec;
}

zoeFeatureVec extractFeature (zoeLabel label) {
	int           i;
	zoeFeatureVec vec = zoeNewFeatureVec();
	
	for (i = 0; i < ANN->vec->size; i++) {
		if (ANN->vec->elem[i]->label == label) {
			zoePushFeatureVec(vec, ANN->vec->elem[i]);
		}
	}
	
	return vec;
}

/*****************************************************************************\
 export functions
\*****************************************************************************/

zoeHash kill_list (const zoeVec genes) {
	zoeCDS        cds1, cds2;
	zoeHash       kill = zoeNewHash();
	int i, j;
		
	/* compare spans using half matrix */
	for (i = 0; i < genes->size; i++) {
		cds1 = genes->elem[i];
		for (j = i + 1; j < genes->size; j++) {
			cds2 = genes->elem[j];
			if (zoeCDSsOverlap(cds1, cds2)) {
				zoeSetHash(kill, cds1->name, (void*)1);
				zoeSetHash(kill, cds2->name, (void*)1);
			}
		}
	}
		
	return kill;
}

void export (void) {
	int         i, j, padding, start, end, leader, trailer;
	zoeCDS     gene /*, check*/;
	zoeVec     genes;
	zoeDNA     dna, tmp;
	zoeFeature exon;
	FILE       * dna_stream;
	FILE       * ann_stream;
	FILE       * tx_stream;
	FILE       * aa_stream;
	/*char         name[256];*/
	
	padding = atoi(zoeOption("-export"));
	if (padding < 3) {
		padding = 3;
		zoeE("setting padding to minimum value of 3\n");
		/* necessary because if stop is not labeled, it may need to be extended */
	}
	
	if ((dna_stream = fopen("export.dna", "w")) == NULL) zoeExit("file open error");
	if ((ann_stream = fopen("export.ann", "w")) == NULL) zoeExit("file open error");
	if ((tx_stream  = fopen("export.tx",  "w")) == NULL) zoeExit("file open error");
	if ((aa_stream  = fopen("export.aa",  "w")) == NULL) zoeExit("file open error");
	
	openData();
	while (getData()) {
		genes = zoeGetGenes(ANN, DNA);	
		
		for (i = 0; i < genes->size; i++) {
			gene = genes->elem[i];
			if (!gene->OK && !zoeOption("-errors-ok")) {
				zoeO("%s %s skipped due to errors\n", DNA->def, gene->name);
				zoeReportCDS(stdout, gene);
				zoeWriteFullCDS(stdout, gene);
				continue;
			}
			
			leader  = gene->start - padding;
			if (leader < 0) leader = 0;
			trailer = gene->end + padding;
			if (trailer > DNA->length) trailer = DNA->length -1;
			
			/* genomic */
			dna = zoeSubseqDNA(gene->name, gene->dna, leader, trailer - leader + 1);
						
			/* annotation coordinate munging */
			for (j = 0; j < gene->exons->size; j++) {
				gene->exons->elem[j]->start -= leader;
				gene->exons->elem[j]->end   -= leader;
			}
			
			/* DNA and annotation strand munging */
			if (gene->strand == '-' && zoeOption("-plus")) {
				for (j = 0; j < gene->exons->size; j++) {
					exon = gene->exons->elem[j];
					
					start = dna->length - exon->end   -1;
					end   = dna->length - exon->start -1;
					
					exon->start  = start;
					exon->end    = end;
					exon->strand = exon->strand == '+' ? '-' : '+';
					
					if (exon->start < 0) zoeExit("aha");
					
				}
				
				tmp = zoeAntiDNA(dna->def, dna);
				zoeDeleteDNA(dna);
				dna = tmp;
			}
			
			/* output */
			zoeWriteDNA(dna_stream, dna);
			fprintf(ann_stream, ">%s\n", dna->def);
			zoeWriteTriteCDS(ann_stream, gene);
			if (gene->tx) zoeWriteDNA(tx_stream, gene->tx);
			else          zoeS(tx_stream, ">%s\n", gene->name);
			if (gene->aa) zoeWriteProtein(aa_stream, gene->aa);
			else          zoeS(aa_stream, ">%s\n", gene->name);
			zoeDeleteCDS(gene);
			zoeDeleteDNA(dna);
		}
		
		/* clean up */
		zoeDeleteVec(genes);
	}
	closeData();
	
}

/*****************************************************************************\
 categorize functions
\*****************************************************************************/

struct GeneRegion {
	zoeVec        regions;
	zoeFeatureVec spans;
};

struct GeneRegion makeRegions(const zoeVec genes, int padding) {
	zoeVec        regions, region;
	zoeFeatureVec spans;
	zoeFeature    span, tx;
	int           i, j, start, end;
	zoeCDS        cds;
	char        * r;
	struct GeneRegion GR;
	
	/* paint the DNA with CDS spans */
	r = zoeCalloc(DNA->length, sizeof(char));
	for (i = 0; i < genes->size; i++) {
		cds = genes->elem[i];
		for (j = cds->start; j <= cds->end; j++) {
			r[j] = 1;
		}
	}
	
	/* create span features */
	spans = zoeNewFeatureVec();
	start = r[0];
	for (i = 1; i < DNA->length; i++) {
		if (r[i]) {
			if (r[i-1] == 0) start = i;
			if (r[i+1] == 0) {
				end = i;
				span = zoeNewFeature(Misc, start, end, '=', 0, 0, 0, 0, NULL/*, NULL*/);
				zoePushFeatureVec(spans, span);
				zoeDeleteFeature(span);
			}
		}
	}
	
	/* find intersection between genes and spans */
	regions = zoeNewVec();
	for (i = 0; i < spans->size; i++) {
		span = spans->elem[i];
		region = zoeNewVec();
		
		for (j = 0; j < genes->size; j++) {
			cds = genes->elem[j];
			
			/*
			if (strcmp(cds->name, "Y48G1C.6") == 0) {
				zoeO("%d..%d %c\n", cds->start, cds->end, cds->strand);
				zoeWriteCDS(stdout, cds);
				zoeWriteFeature(stdout, span);
				zoeExit("got it");
			}
			*/
			
			tx  = zoeNewFeature(Coding, cds->start, cds->end, '=', 0, 0, 0, 0, NULL/*, NULL*/);
			
			if (zoeFeaturesOverlap(span, tx)) zoePushVec(region, cds);
			
			zoeDeleteFeature(tx);
		}
		
		zoePushVec(regions, region);
	}
	
	GR.regions = regions;
	GR.spans   = spans;
	return GR;
}

char categorize_locus (const char * name, zoeVec genes) {
	int i, j, errors, warnings, share, overlap;
	zoeCDS gene1, gene2;
	
	/* categories:               char
		unique transcript         u
		alternate transcripts     a
		overlapping transcripts   o
		errors                    e
		warnings                  w
	*/
	
	warnings = 0;
	errors   = 0;
	share    = 0;
	overlap  = 0;
	
	for (i = 0; i < genes->size; i++) {
		gene1 = genes->elem[i];
		if (gene1->errors->size)   errors++;
		if (gene1->warnings->size) warnings++;
		
		for (j = i+1; j < genes->size; j++) {
			gene2 = genes->elem[j];
			if (zoeCDSsOverlap(gene1, gene2))       overlap++;
			if (zoeCDSsShareSequence(gene1, gene2)) share++;
		}
	}
	
	if (errors)            return 'e';
	if (warnings)          return 'w';
	if (overlap && share)  return 'a';
	if (overlap && !share) return 'o';
	if (genes->size == 1)  return 'u';
	if (genes->size > 1)   return 'o'; /* 1 bp apart?, put it here for now */
	
	zoeE("genes=%d, overlap=%d, share=%d\n", genes->size, overlap, share);
	for (i = 0; i < genes->size; i++) {
		gene1 = genes->elem[i];
		zoeE("%s\n", gene1->name);
	}
	zoeExit("error in categorize_locus");
	
	return 'X';
}


void categorize (void) {
	int               i, j, k, padding, d, start, end, span_count = 0;
	struct GeneRegion GR;
	zoeVec            genes, loci;
	zoeFeature        span, prev_span, post_span, exon;
	zoeFeatureVec     exons;
	zoeFeatureTable   ann;
	zoeDNA            dna;
	zoeCDS            cds;
	char              category, name[1024];
	FILE            * e_dna, * e_ann;
	FILE            * w_dna, * w_ann;
	FILE            * a_dna, * a_ann;
	FILE            * o_dna, * o_ann;
	FILE            * u_dna, * u_ann;
	FILE            * dna_stream, * ann_stream;
	
	padding = atoi(zoeOption("-categorize"));

	if ((e_dna = fopen("err.dna", "w")) == NULL) zoeExit("file open error");
	if ((e_ann = fopen("err.ann", "w")) == NULL) zoeExit("file open error");
	if ((w_dna = fopen("wrn.dna", "w")) == NULL) zoeExit("file open error");
	if ((w_ann = fopen("wrn.ann", "w")) == NULL) zoeExit("file open error");
	if ((a_dna = fopen("alt.dna", "w")) == NULL) zoeExit("file open error");
	if ((a_ann = fopen("alt.ann", "w")) == NULL) zoeExit("file open error");
	if ((o_dna = fopen("olp.dna", "w")) == NULL) zoeExit("file open error");
	if ((o_ann = fopen("olp.ann", "w")) == NULL) zoeExit("file open error");
	if ((u_dna = fopen("uni.dna", "w")) == NULL) zoeExit("file open error");
	if ((u_ann = fopen("uni.ann", "w")) == NULL) zoeExit("file open error");
	
	openData();
	while (getData()) {
		genes = zoeGetGenes(ANN, DNA);
		GR = makeRegions(genes, padding);
		
		for (i = 0; i < GR.spans->size; i++) {
			span = GR.spans->elem[i];
			loci = GR.regions->elem[i];
			
			/* extend span start */
			if (i > 0) {
				prev_span = GR.spans->elem[i-1];
				d = (float)(span->start - prev_span->end + 1) / 2;
				if (d > padding) d = padding;
				start = span->start - d;
			} else {
				start = span->start - padding;
				if (start < 1) start = 0;
			}
						
			/* extend span end */
			if (i < GR.spans->size -1) {
				post_span = GR.spans->elem[i+1];
				d = (float)(post_span->start - span->end + 1) / 2;
				if (d > padding) d = padding;
				end = span->end + d;
			} else {
				end = span->end + padding;
				if (end >= DNA->length) end = DNA->length -1;
			}
			
			if (start < 0 || end < 0) zoeExit("negatory");
			
			span_count++;
			sprintf(name, "region-%d (%s %d..%d)", span_count, DNA->def, start, end);
			
			/* categorize */
			category = categorize_locus(name, loci);			
			switch (category) {
				case 'e': dna_stream = e_dna; ann_stream = e_ann; break;
				case 'w': dna_stream = w_dna; ann_stream = w_ann; break;
				case 'a': dna_stream = a_dna; ann_stream = a_ann; break;
				case 'o': dna_stream = o_dna; ann_stream = o_ann; break;
				case 'u': dna_stream = u_dna; ann_stream = u_ann; break;
				default:  dna_stream = NULL;  ann_stream = NULL;
			}
						
			/* output DNA */
			dna = zoeSubseqDNA(name, DNA, start, end - start + 1);
			zoeWriteDNA(dna_stream, dna);
			zoeDeleteDNA(dna);
			
			/* output ANN */
			exons = zoeNewFeatureVec();
			for (j = 0; j < loci->size; j++) {
				cds = loci->elem[j];
				for (k = 0; k < cds->exons->size; k++) {
					exon = cds->exons->elem[k];
					exon->start -= start;
					exon->end   -= start;
					if (!zoeVerifyFeature(exon)) {
						zoeWarn("adjusting coordinates in categorize killed exon");
						zoeExit("%s", cds->name);
					}
					zoePushFeatureVec(exons, exon);
				}
			}
			
			ann = zoeNewFeatureTable(name, exons);
			zoeWriteTriteFeatureTable(ann_stream, ann);
			
			zoeDeleteFeatureVec(exons);
			zoeDeleteFeatureTable(ann);
		}
		
		/* clean up */
		zoeDeleteVec(GR.regions);
		zoeDeleteFeatureVec(GR.spans);
		for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
		zoeDeleteVec(genes);
	}
	closeData();
}

void export_feature (void) {
	int         i, start, end, upstream, downstream;
	zoeDNA     dna;
	zoeFeature f;
	zoeLabel   label = zoeText2Label(zoeOption("-export-feature"));
	FILE       * dna_stream;
	FILE       * ann_stream;
	
	if (zoeOption("-upstream")) upstream = atoi(zoeOption("-upstream"));
	else                        upstream = 0;
	if (zoeOption("-downstream")) downstream = atoi(zoeOption("-downstream"));
	else                          downstream = 0;
	
	if ((dna_stream = fopen("export-feature.dna", "w")) == NULL) zoeExit("file open error");
	if ((ann_stream = fopen("export-feature.ann", "w")) == NULL) zoeExit("file open error");

	
	openData();
	while (getData()) {
		for (i = 0; i < ANN->vec->size; i++) {
			f = ANN->vec->elem[i];
			if (f->label != label) continue;
			
			start = f->start - upstream;
			end   = f->end + downstream;
			if (start < 0 || end >= DNA->length) {
				zoeWarn("feature out of bounds, skipping");
				zoeWriteFeature(stderr, f);
				continue;
			}
			dna = zoeSubseqDNA(f->group, DNA, start, end - start +1);
			zoeWriteDNA(dna_stream, dna);
			
			
			zoeWriteFeature(ann_stream, f);
			zoeDeleteDNA(dna);
		}
	}
	closeData();
	
}

/*****************************************************************************\
 score genes
\*****************************************************************************/

void score_genes (void) {
	int           i, errors_ok;
	zoeCDS        cds;
	zoeVec        genes;
	zoeHMM        hmm;
	zoeTrellis    t;
	
	hmm = zoeGetHMM(zoeOption("-score-genes"));
	zoeSetTrellisMeter(0);
	
	if (zoeOption("-errors-ok")) errors_ok = 1;
	else                         errors_ok = 0;
	
	openData();
	while (getData()) {
		t = zoeNewTrellis(DNA, hmm, NULL);
	
		genes = zoeGetGenes(ANN, DNA);
		printf(">%s\n", DNA->def);
		for (i = 0; i < genes->size; i++) {
			cds = genes->elem[i];
			zoeScoreCDS(t, cds, 0, errors_ok);
			zoeWriteFullCDS(stdout, cds);
		}
		
		for (i = 0; i < genes->size; i++) zoeDeleteCDS((zoeCDS)genes->elem[i]);
		zoeDeleteVec(genes);
		zoeDeleteTrellis(t);
	}
	closeData();
}


/*****************************************************************************\
 filter genes
\*****************************************************************************/

void filter_genes (void) {
	int           i;
	zoeCDS        cds;
	zoeVec        genes;
	zoeHMM        hmm;
	zoeTrellis    t;
	float         min_score;
	int           min_length;
	
	min_score  = atof(zoeOption("-min-score"));
	min_length = atoi(zoeOption("-min-length"));
	
	hmm = zoeGetHMM(zoeOption("-filter-genes"));
	zoeSetTrellisMeter(0);
	
	openData();
	while (getData()) {
		t = zoeNewTrellis(DNA, hmm, NULL);
	
		genes = zoeGetGenes(ANN, DNA);
		printf(">%s\n", DNA->def);
		for (i = 0; i < genes->size; i++) {
			cds = genes->elem[i];
			zoeScoreCDS(t, cds, 0, 0);
			/*printf("%.1f\t%d\n", cds->score, cds->tx->length);*/
			if (cds->score > min_score && cds->tx->length > min_length) {
				zoeWriteFullCDS(stdout, cds);
			}
		}
		
		for (i = 0; i < genes->size; i++) zoeDeleteCDS((zoeCDS)genes->elem[i]);
		zoeDeleteVec(genes);
		zoeDeleteTrellis(t);
	}
	closeData();
}

/*****************************************************************************\
 compare genes
\*****************************************************************************/

struct nt_stats {
	int a_only;
	int b_only;
	int shared;
};

struct nt_stats calculate_nt_stats (int length, zoeFeatureVec v1, zoeFeatureVec v2) {
	int             i, j, *a, *b;
	zoeFeature      f;
	struct nt_stats stats;
		
	stats.a_only = 0;
	stats.b_only = 0;
	stats.shared = 0;
	
	a = zoeCalloc(length * sizeof(int), 1);
	b = zoeCalloc(length * sizeof(int), 1);
	
	for (i = 0; i < v1->size; i++) {
		f = v1->elem[i];
		for (j = f->start; j <= f->end; j++) a[j] = 1;
	}
	
	for (i = 0; i < v2->size; i++) {
		f = v2->elem[i];
		for (j = f->start; j <= f->end; j++) b[j] = 1;
	}
	
	for (i = 0; i < length; i++) {
		if (a[i] && b[i]) stats.shared++;
		else if (a[i])    stats.a_only++;
		else if (b[i])    stats.b_only++;
	}
	
	zoeFree(a);
	zoeFree(b);
	
	return stats;
}

struct f_stats {
	int shared[zoeLABELS];
	int a_only[zoeLABELS];
	int b_only[zoeLABELS];
};

struct f_stats calculate_feature_stats (int length, zoeFeatureVec v1, zoeFeatureVec v2) {
	int               i, *a, *b;
	zoeFeature        f;
	struct f_stats    stats;
	
	for (i = 0; i < zoeLABELS; i++) {
		stats.shared[i] = 0;
		stats.a_only[i] = 0;
		stats.b_only[i] = 0;
	}
	
	a = zoeCalloc(length * sizeof(zoeLabel), 1);
	b = zoeCalloc(length * sizeof(zoeLabel), 1);
	
	for (i = 0; i < v1->size; i++) {
		f = v1->elem[i];
		switch (f->label) {
			case Einit: a[f->start] = Start;    a[f->end] = Donor; break;
			case Eterm: a[f->start] = Acceptor; a[f->end] = Stop;  break;
			case Exon:  a[f->start] = Acceptor; a[f->end] = Donor; break;
			case Esngl: a[f->start] = Start;    a[f->end] = Stop;  break;
			default: break;
		}
	}
	
	for (i = 0; i < v2->size; i++) {
		f = v2->elem[i];
		switch (f->label) {
			case Einit: b[f->start] = Start;    b[f->end] = Donor; break;
			case Eterm: b[f->start] = Acceptor; b[f->end] = Stop;  break;
			case Exon:  b[f->start] = Acceptor; b[f->end] = Donor; break;
			case Esngl: b[f->start] = Start;    b[f->end] = Stop;  break;
			default: break;
		}
	}
	
	for (i = 0; i < length; i++) {
		if (a[i] == 0 && b[i] == 0) continue;
		else if (a[i] == b[i])  stats.shared[a[i]]++;
		else if (a[i] && !b[i]) stats.a_only[a[i]]++;
		else if (b[i] && !a[i]) stats.b_only[b[i]]++;
		else {
			stats.a_only[a[i]]++;
			stats.b_only[b[i]]++;
		}
	}
	
	zoeFree(a);
	zoeFree(b);
	return stats;
}

struct f_stats calculate_exon_stats (zoeFeatureVec v1, zoeFeatureVec v2) {
	int               i;
	zoeFeature        f;
	struct f_stats    stats;
	char              string[64];
	zoeHash           h1, h2;
	zoeTVec           keys;
	char            * key;
	void            * val;
	
	for (i = 0; i < zoeLABELS; i++) {
		stats.shared[i] = 0;
		stats.a_only[i] = 0;
		stats.b_only[i] = 0;
	}
	
	h1 = zoeNewHash();
	for (i = 0; i < v1->size; i++) {
		f = v1->elem[i];
		sprintf(string, "%d %d %c", f->start, f->end, f->strand);
		zoeSetHash(h1, string, (void*)1);
	}
	
	h2 = zoeNewHash();
	for (i = 0; i < v2->size; i++) {
		f = v2->elem[i];
		sprintf(string, "%d %d %c", f->start, f->end, f->strand);
		zoeSetHash(h2, string, (void*)1);
	}
	
	keys = zoeKeysOfHash(h1);
	for (i = 0; i < keys->size; i++) {
		key = keys->elem[i];
		val = zoeGetHash(h2, key);
		if (val) {
			stats.shared[(size_t)Exon]++;
		} else {
			stats.a_only[(size_t)Exon]++;
		}
	}
	zoeDeleteTVec(keys);
	
	keys = zoeKeysOfHash(h2);
	for (i = 0; i < keys->size; i++) {
		key = keys->elem[i];
		val = zoeGetHash(h1, key);
		if (val) {
		} else {
			stats.b_only[(size_t)Exon]++;
		}
	}
	zoeDeleteTVec(keys);
	
	zoeDeleteHash(h1);
	zoeDeleteHash(h2);
		
	return stats;
}

struct g_stats {
	int shared;
	int a_only;
	int b_only;
};

struct g_stats calculate_gene_stats (zoeVec v1, zoeVec v2) {
	int            i, j, found;
	zoeCDS         gene1, gene2;
	struct g_stats stats;
	
	stats.shared = 0;
	stats.a_only = 0;
	stats.b_only = 0;
	
	for (i = 0; i < v1->size; i++) {
		gene1 = v1->elem[i];
		found = 0;
		for (j = 0; j < v2->size; j++) {
			gene2 = v2->elem[j];
			if (zoeCDScmp(gene1, gene2) == 0) {
				stats.shared++;
				found = 1;
				break;
			}
		}
		if (!found) stats.a_only++;
	}
	
	for (j = 0; j < v2->size; j++) {
		gene2 = v2->elem[j];
		found = 0;
		for (i = 0; i < v1->size; i++) {
			gene1 = v1->elem[i];
			if (zoeCDScmp(gene2, gene1) == 0) {
				found = 1;
				break;
			}
		}
		if (!found) stats.b_only++;
	}
		
	return stats;
}

static float sensitivity (int shared, int a_only, int b_only) {
	int num, den;
	
	num = shared;
	den = shared + a_only;
	
	if (den == 0) return 0;
	else          return (float)num / (float)den;
}

static float specificity (int shared, int a_oly, int b_only) {
	int num, den;

	num = shared;
	den = shared + b_only;
	
	if (den == 0) return 0;
	else          return (float)num / (float)den;
}

void compare_genes (void) {
	int             i, j;
	zoeVec          genes1, genes2;
	zoeCDS          gene;
	zoeFeatureVec   exons1, exons2;
	zoeFeature      exon;
	struct nt_stats nt, NT;
	struct f_stats  fs, FS, es, ES;
	struct g_stats  gs, GS;
	char            name[16];

	/* init */
	NT.a_only = 0; NT.b_only = 0; NT.shared = 0;
	GS.a_only = 0; GS.b_only = 0; GS.shared = 0;
	for (i = 0; i < zoeLABELS; i++) {
		FS.shared[i] = 0; FS.a_only[i] = 0; FS.b_only[i] = 0;
		ES.shared[i] = 0; ES.a_only[i] = 0; ES.b_only[i] = 0;
	}
	
	/* processing loop */
	openData();
	while (getData()) {
		
		genes1 = zoeGetGenes(ANN, DNA);
		genes2 = zoeGetGenes(PRE, DNA);
				
		/* collect all exons */
		exons1 = zoeNewFeatureVec();
		exons2 = zoeNewFeatureVec();
		for (i = 0; i < genes1->size; i++) {
			gene = genes1->elem[i];
			for (j = 0; j < gene->exons->size; j++) {
				exon = gene->exons->elem[j];
				zoePushFeatureVec(exons1, exon);
			}
		}
		for (i = 0; i < genes2->size; i++) {
			gene = genes2->elem[i];
			for (j = 0; j < gene->exons->size; j++) {
				exon = gene->exons->elem[j];
				zoePushFeatureVec(exons2, exon);
			}
		}
		
		/* nucleotide stats */
		nt = calculate_nt_stats(DNA->length, exons1, exons2);
		NT.a_only += nt.a_only;
		NT.b_only += nt.b_only;
		NT.shared += nt.shared;
		
		/* feature stats */
		fs = calculate_feature_stats(DNA->length, exons1, exons2);
		for (i = 0; i < zoeLABELS; i++) {
			FS.shared[i] += fs.shared[i];
			FS.a_only[i] += fs.a_only[i];
			FS.b_only[i] += fs.b_only[i];
		}
		
		/* exon stats - exact Einit, etc */
		es = calculate_exon_stats(exons1, exons2);
		for (i = 0; i < zoeLABELS; i++) {
			ES.shared[i] += es.shared[i];
			ES.a_only[i] += es.a_only[i];
			ES.b_only[i] += es.b_only[i];
		}
		
		/* gene stats */
		gs = calculate_gene_stats(genes1, genes2);
		GS.shared += gs.shared;
		GS.a_only += gs.a_only;
		GS.b_only += gs.b_only;
		
		/* protein & tx stats - would be neat to do N-W here */
		
		/* verbose, sequence reporting */
		if (zoeOption("-details")) {
			zoeO("%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
				DNA->def, 
				(int) (100 * sensitivity(nt.shared, nt.a_only, nt.b_only)),
				(int) (100 * specificity(nt.shared, nt.a_only, nt.b_only)),

				fs.shared[Acceptor],
				fs.a_only[Acceptor],
				fs.b_only[Acceptor],
				
				fs.shared[Donor],
				fs.a_only[Donor],
				fs.b_only[Donor],
				
				fs.shared[Start],
				fs.a_only[Start],
				fs.b_only[Start],
				
				fs.shared[Stop],
				fs.a_only[Stop],
				fs.b_only[Stop]
			);
		}
		
		/* clean up */
		for (i = 0; i < genes1->size; i++) zoeDeleteCDS(genes1->elem[i]);
		zoeDeleteVec(genes1);
		for (i = 0; i < genes2->size; i++) zoeDeleteCDS(genes2->elem[i]);
		zoeDeleteVec(genes2);
		zoeDeleteFeatureVec(exons1);
		zoeDeleteFeatureVec(exons2);
	}
	closeData();
	
	/* Global Ouput */
	if (zoeOption("-details")) return;
	zoeO("Nucleotide: %.3f %.3f\n",
		sensitivity(NT.shared, NT.a_only, NT.b_only),
		specificity(NT.shared, NT.a_only, NT.b_only));
	zoeO("Gene: %.3f %.3f\n",
		sensitivity(GS.shared, GS.a_only, GS.b_only),
		specificity(GS.shared, GS.a_only, GS.b_only));
	for (i = 0; i < zoeLABELS; i++) {
		if (FS.shared[i] || FS.a_only[i] || FS.b_only[i]) {
			zoeLabel2Text(i, name);
			zoeO("%s: %.3f %.3f\n", name,
				sensitivity(FS.shared[i], FS.a_only[i], FS.b_only[i]),
				specificity(FS.shared[i], FS.a_only[i], FS.b_only[i]));
		}
	}
	for (i = 0; i < zoeLABELS; i++) {
		if (ES.shared[i] || ES.a_only[i] || ES.b_only[i]) {
			zoeLabel2Text(i, name);
			zoeO("%s: %.3f %.3f\n", name,
				sensitivity(ES.shared[i], ES.a_only[i], ES.b_only[i]),
				specificity(ES.shared[i], ES.a_only[i], ES.b_only[i]));
		}
	}
	
}

/*****************************************************************************\
 ace format
\*****************************************************************************/

void ace_format (void) {
	int           i, j;
	coor_t        start, end;
	zoeFeature    exon;
	zoeVec        genes;
	zoeCDS        cds;
	zoeFeatureVec exons;
	char          * gene_method = zoeOption("-gene-method");
	char          * extra = zoeOption("-extra");

	if (!zoeOption("-gene-method")) zoeExit("-gene-method required");
	if (extra == NULL) extra = "";
	
	openData();
	while (getData()) {
		if (zoeOption("-dna")) {
			zoeO("\nDNA %s\n", DNA->def);
			zoeO("%s\n", DNA->seq);
		
			zoeO("\nSequence %s\n", DNA->def);
			zoeO("Genomic_canonical\n");
		}
		
		genes = zoeGetGenes(ANN, DNA);
		exons = zoeNewFeatureVec();
		for (i = 0; i < genes->size; i++) {
			cds = genes->elem[i];
			
			if (cds->strand == '+') {
				start = cds->start +1;
				end   = cds->end   +1;
			} else {
				start = cds->end   +1;
				end   = cds->start +1;
			}
			
			zoeO("\nSequence %s\n", DNA->def);
			zoeO("Subsequence %s%s %d %d\n", cds->name, extra, start, end);
			zoeO("\nSequence %s%s\n", cds->name, extra);
			
			if (cds->strand == '+') {
				for (j = 0; j < cds->exons->size; j++) {
					exon = cds->exons->elem[j];
					/*zoePushFeatureVec(exons, exon);*/
					start = exon->start - cds->start +1;
					end   = exon->end   - cds->start +1;
					zoeO("Source_Exons %d %d\n", start, end);
				}
			} else {
				for (j = 0; j < cds->exons->size; j++) {
					exon = cds->exons->elem[j];
					/*zoePushFeatureVec(exons, exon);*/
					end   = cds->end - exon->start +1;
					start = cds->end - exon->end   +1;
					zoeO("Source_Exons %d %d\n", start, end);
				}
			}

			zoeO("CDS\nMethod %s\n", gene_method);
			zoeDeleteCDS(cds);
		}
		
		zoeDeleteVec(genes);
		zoeDeleteFeatureVec(exons);
	}
	closeData();
}

/*****************************************************************************\
 generic data processing functions
\*****************************************************************************/

void openData (void) {
	DNA_stream = zoeOpenFile(DNA_file);
	ANN_stream = zoeOpenFile(ANN_file);
	if (PREDICTIONS) PRE_stream = zoeOpenFile(PRE_file);
}

void closeData (void) {

	zoeCloseFile(DNA_stream);
	zoeCloseFile(ANN_stream);
	if (PREDICTIONS) zoeCloseFile(PRE_stream);
	
	if (DNA) {
		zoeDeleteDNA(DNA);
		DNA = NULL;
	}
	
	if (ANTI) {
		zoeDeleteDNA(ANTI);
		ANTI = NULL;
	}
	
	if (ANN) {
		zoeDeleteFeatureTable(ANN);
		ANN = NULL;
	}
	
	if (PRE) {
		zoeDeleteFeatureTable(PRE);
		PRE = NULL;
	}
	
}

int getData (void) {
	zoeFastaFile ff;

	/* read new dna and ann */
	if ((ff = zoeReadFastaFile(DNA_stream.stream)) == NULL) return 0;
	
	/* delete old dna and ann objects */
	if (DNA)  zoeDeleteDNA(DNA);
	if (ANTI) zoeDeleteDNA(ANTI);
	if (ANN)  zoeDeleteFeatureTable(ANN);
	if (PRE)  zoeDeleteFeatureTable(PRE);
	
	DNA = zoeNewDNA(ff->def, ff->seq);
	zoeDeleteFastaFile(ff);
	if (ANTI_REQUIRED) ANTI = zoeAntiDNA(DNA->def, DNA);
	ANN = zoeReadFeatureTable(ANN_stream.stream);
	if (PREDICTIONS) PRE = zoeReadFeatureTable(PRE_stream.stream);
	
	/* lowercase */
	if (zoeOption("-lcfilter")) {
		zoeLCfilter(DNA);
		if (ANTI_REQUIRED) zoeLCfilter(ANTI);
	}
	if (zoeOption("-lcsmooth")) {
		zoeLCsmooth(DNA, 10, 10, 100);
		if (ANTI_REQUIRED) zoeLCsmooth(ANTI, 10, 10, 100);
	}
	
	return 1;
}

void help (void) {

zoeM(stdout, 47,

"\nThe general form of the fathom command line is:\n",
"    fathom <ZFF file> <FASTA file> <command> [options]\n",
"There are many commands, and not all are listed in the usage statement.",
"The complete list is given farther below.",
"\nZFF file:\n",
"    The annotation file must be formatted in ZFF, which is an odd mixture",
"    of FASTA and GFF.",
"\nFASTA file:\n",
"    The DNA must be in FASTA format. The definition line of each ZFF and",
"    FASTA file must be identical and in the same order. If the sequences",
"    have been masked with lowercase letters, include -lcmask.",
"\nCommands:\n",
"-validate       Checks the annotation for warnings and errors and reports",
"                a short summary to stdout.",
"-gene-stats     Reports various statistics about the dataset.",
"-export         Creates 4 files export.{ann,dna,aa,tx}. The coordinates and",
"                dna file depend on the value of <int>. The -plus option",
"                converts the sequence to plus strand.",
"-categorize     Categorizes genomic regions into those that contain errors",
"                (err), warnings (wrn), alternate forms (alt), overlapping",
"                genes (olp), and unique genes (uni). Typically, only the",
"                unique genes are used for training and testing. The value",
"                of <i> limits the intergenic sequence at the ends.",
"-split          Splits the dataset into several pieces (-number) by fraction",
"                (-training) or by GC. For -training, the random seed may be",
"                set via -randseed.",
"-ace-format     Converts annotation to ACEDB format. A display method must",
"                be specified by -gene-method. To convert the dna to ACEDB",
"                format, include the -dna option.",
"-compare-genes  Calculates sensitivity and specificity at the nucleotide",
"                exon, and gene levels, as well as various sequence models.",
"-extract        Extracts sequence around a specific sequence feature. The",
"                -offset and -length determine the region.",
"-reformat       Converts trite annotation to full annotation.",
"-anti           Reverse-complements the sequence and annotation.",
"-fuse <int>     Merges sequences and annotation.",
"-lcfilter       Converts lowercase to N.",
"-lcsmooth       Smooths small islands of N or non-N.",
"-min-intron     Sets the warning value for minimum intron length (default 30).",
"-max-intron     (default 100000)",
"-min-exon       (default 6)",
"-max-exon       (default 50000)",
"-min-gene       (default 150)",
"-max-gene       (default 200000)",
"-min-cds        (default 150)",
"-max-cds        (default 50000)",
"\nTry 'fathom -help | more' if the output has scrolled off your screen"
);

exit(0);

}
