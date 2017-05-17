/*****************************************************************************\
 snap.c

Semi-HMM-based Nucleic Acid Parser, a gene prediction program

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

int    file_is_isochore (const char *);
float  gc_fraction (const zoeDNA);
int    dna_is_ok (const zoeDNA);
zoeVec parse_dna (const zoeHMM, const zoeDNA, zoeFeatureTable);

void   ace_output (const zoeDNA, const zoeVec);
void   gff_output (const zoeDNA, const zoeVec);
void   zoe_output (const zoeDNA, const zoeVec);
void   help (void);

score_t SNAP_OVERLAP   = 200;
coor_t  SNAP_MIN_CDS   = 180;
score_t SNAP_MIN_SCORE = -1000;
char * ZOE = NULL; /* environment variable */


static char usage[]  = "\n\
SNAP - Semi-HMM-based Nucleic Acid Parser (version 2006-07-28)\n\n\
usage: snap [options] <HMM file> <FASTA file> [options]\n\
options:\n\
  -help           report useful information\n\
  -lcmask         treat lowercase as N\n\
  -plus           predict on plus strand only\n\
  -minus          predict on minus strand only\n\
  -gff            output annotation as GFF\n\
  -ace            output annotation as ACED\n\
  -quiet          do not send progress to STDERR\n\
  -aa <file>      create FASTA file of proteins\n\
  -tx <file>      create FASTA file of transcripts\n\
  -xdef <file>    external definitions\n\
  -name <string>  name for the gene [default snap]\n\
";

/*
static char help[]
*/

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char *argv[]) {
	zoeFile         dna_file;
	zoeFastaFile    fasta;
	zoeDNA          dna;
	zoeFeatureTable xdef = NULL;
	zoeHMM          hmm;
	zoeIsochore     iso;
	zoeCDS          gene;
	zoeVec          genes;
	int             label, i;
	char            option[16], name[16];
	FILE          * aa_stream = NULL;
	FILE          * tx_stream = NULL;
	FILE          * xd_stream  = NULL;
	
	/* set the program name */
	zoeSetProgramName(argv[0]);
	
	/* general options */
	zoeSetOption("-help",    0);
	zoeSetOption("-lcmask",  0);
	zoeSetOption("-plus",    0);
	zoeSetOption("-minus",   0);
	zoeSetOption("-quiet",   0);
	zoeSetOption("-gff",     0);
	zoeSetOption("-ace",     0);
	zoeSetOption("-aa",      1);
	zoeSetOption("-tx",      1);
	zoeSetOption("-xdef",    1);
	
	/* unadvertised options for my own use/testing */
	zoeSetOption("-name",      1);
	zoeSetOption("-overlap",   1);
	zoeSetOption("-min-cds",   1);
	zoeSetOption("-min-score", 1);
	zoeSetOption("-flatN",     0);
	zoeSetOption("-boostN",    0);
	zoeSetOption("-debug",     0);
	zoeSetOption("-xdebug",    0);
	zoeSetOption("-info",      0);
	zoeSetOption("-extra",     0);

	
	/* Nscore parameters - what to do if there is an N */
	for (label = 0; label < zoeLABELS; label++) {
		zoeLabel2Text(label, name);
		sprintf(option, "-N%s", name);
		zoeSetOption(option, 1);
	}
	
	/* Ascore parameters - what to modify all scores by */
	for (label = 0; label < zoeLABELS; label++) {
		zoeLabel2Text(label, name);
		sprintf(option, "-A%s", name);
		zoeSetOption(option, 1);
	}
	
	zoeParseOptions(&argc, argv);
	if (zoeOption("-help")) help();
		
	/* usage */
	if (argc != 3) {
		zoeE("%s", usage);
		exit(1);
	}
	
	/* ZOE environment variable */
	ZOE = getenv("ZOE");
	
	/* strand */
	if (zoeOption("-plus") && zoeOption("-minus")) {
		zoeExit("-plus or -minus. Omit for both strands.");
	}
	
	/* quiet */
	if (zoeOption("-quiet")) zoeSetTrellisMeter(0);
	
	/* others */
	if (zoeOption("-overlap")) SNAP_OVERLAP = atof(zoeOption("-overlap"));
	if (zoeOption("-min-cds")) SNAP_MIN_CDS = atoi(zoeOption("-min-cds"));
	if (zoeOption("-min-score")) SNAP_MIN_SCORE = atof(zoeOption("-min-score"));
	
	/* -flatN and -boostN */
	if (zoeOption("-flatN") && zoeOption("-boostN")) {
		zoeExit("-flatN and -boostN are mutually incompatible");
	} else if (zoeOption("-flatN")) {
		for (label = 0; label < zoeLABELS; label++) zoeSetNscore(label, -1);
	} else if (zoeOption("-boostN")) {
		zoeSetNscore(Acceptor, -10);
		zoeSetNscore(Donor,    -10);
		zoeSetNscore(Start,    -10);
		zoeSetNscore(Stop,     -10);
		zoeSetNscore(Coding,   -10);
		zoeSetNscore(Repeat,   +10);
		zoeSetNscore(Inter,    +10);
		zoeSetNscore(Intron,   +10);
	}
	
	/* set individual Nscore values */
	for (label = 0; label < zoeLABELS; label++) {
		zoeLabel2Text(label, name);
		sprintf(option, "-N%s", name);
		if (zoeOption(option)) zoeSetNscore(label, atof(zoeOption(option)));
	}
	
	/* set individual Ascore values */
	for (label = 0; label < zoeLABELS; label++) {
		zoeLabel2Text(label, name);
		sprintf(option, "-A%s", name);
		if (zoeOption(option)) zoeSetAscore(label, atof(zoeOption(option)));
	}
	
	/* HMM and isochore */
	if (file_is_isochore(argv[1])) {
		iso = zoeGetIsochore(argv[1]);
		hmm = NULL;
	} else {
		iso = NULL;
		hmm = zoeGetHMM(argv[1]);
	}
	
	/* more Nscore stuff */
	for (label = 0; label < zoeLABELS; label++) {
		zoeLabel2Text(label, name);
		sprintf(option, "-%s", name);
		if (zoeOption(option) && hmm->mmap[label] == NULL) {
			zoeExit("You have set %s, but the HMM (%s) has no %s model\n",
				option, argv[1], name);
		}
	}
	
	/* aa and tx */
	if (zoeOption("-aa")) {
		if ((aa_stream = fopen(zoeOption("-aa"), "w")) == NULL) {
			zoeExit("error opening -aa file");
		}
	}
	if (zoeOption("-tx")) {
		if ((tx_stream = fopen(zoeOption("-tx"), "w")) == NULL) {
			zoeExit("error opening -tx file");
		}
	}
	
	/* Fasta */
	dna_file = zoeOpenFile(argv[2]);

	/* Xdef */
	if (zoeOption("-xdef")) {
		if ((xd_stream = fopen(zoeOption("-xdef"), "r")) == NULL)
			zoeExit("error opening xdef file %s", zoeOption("-xdef"));
	}
	
	/***************\
		Main Loop
	\***************/
	while ((fasta = zoeReadFastaFile(dna_file.stream)) != NULL) {
	
		/* DNA */
		dna = zoeNewDNA(fasta->def, fasta->seq);
		zoeDeleteFastaFile(fasta);
		if (zoeOption("-lcmask")) {
			zoeLCsmooth(dna, 10, 10,100);
		}
				
		/* Xdef */
		if (zoeOption("-xdef")) {
			xdef = zoeReadFeatureTable(xd_stream);
		}
		
		/* must set isochore hmm if isochores in use */
		if (iso) {
			hmm = zoeSelectIsochore(iso, gc_fraction(dna));
		}
		
		/* check DNA */
		if (dna_is_ok(dna)) genes = parse_dna(hmm, dna, xdef);
		else                genes = zoeNewVec();
		
		/* annotation output */
		if      (zoeOption("-gff")) gff_output(dna, genes);
		else if (zoeOption("-ace")) ace_output(dna, genes);
		else                        zoe_output(dna, genes);
		
		/* sequence output */
		for (i = 0; i < genes->size; i++) {
			gene = genes->elem[i];
			if (zoeOption("-aa")) zoeWriteProtein(aa_stream, gene->aa);
			if (zoeOption("-tx")) zoeWriteDNA(tx_stream, gene->tx);
		}

		/* clean up */
		for (i = 0; i < genes->size; i++) zoeDeleteCDS(genes->elem[i]);
		zoeDeleteVec(genes);
		if (zoeOption("-xdef")) zoeDeleteFeatureTable(xdef);
		zoeDeleteDNA(dna);
	}
	
	if (aa_stream) fclose(aa_stream);
	if (tx_stream) fclose(tx_stream);
	if (xd_stream) fclose(xd_stream);
	zoeCloseFile(dna_file);
	
	if (iso) zoeDeleteIsochore(iso);
	else     zoeDeleteHMM(hmm);
		
	return 0;
}

int cmp_genes_by_score (const zoeCDS a, const zoeCDS b) {
	return b->score - a->score;
}

int cmp_genes_ptr_by_score (const void * a, const void * b) {
	return cmp_genes_by_score( *(zoeCDS *)a, *(zoeCDS *)b );
}

int cmp_genes_by_start (const zoeCDS a, const zoeCDS b) {
	return a->start - b->start;
}

int cmp_genes_ptr_by_start (const void * a, const void * b) {
	return cmp_genes_by_start( *(zoeCDS *)a, *(zoeCDS *)b );
}

void edit_names (zoeFeatureVec vec, const char * name) {
	int        i;
	zoeFeature f;
	
	for (i = 0; i < vec->size; i++) { 
		f = vec->elem[i];
		if (f->group) zoeFree(f->group);
		f->group = zoeMalloc(strlen(name) +1);
		strcpy(f->group, name);
	}
}

int file_is_isochore (const char * name) {
	FILE * file;
	char   type[64];
	char   path[1024];
	
	file = fopen(name, "r");
	if (file == NULL) {
		sprintf(path, "%s/HMM/%s", ZOE, name);
		file = fopen(path, "r");
		if (file == NULL) zoeExit("error opening file (%s)", path);
	}
	
	if (fscanf(file, "%s", type) != 1) zoeExit("error checking HMM file");
	fclose(file);
	
	if (strcmp(type, "zoeHMM") == 0) return 0;
	if (strcmp(type, "zoeIsochore") == 0) return 1;
	
	zoeExit("unrecognized parameter file");
	
	return 0;
}

float gc_fraction (const zoeDNA dna) {
	int i, c[5];
	
	for (i = 0; i < 5; i++) c[i] = 0;
	for (i = 0; i < dna->length; i++) c[(int)dna->s5[i]]++;
	
	return (float)(c[1]+c[2]) / (float)(c[1]+c[2]+c[3]+c[4]);
}

int dna_is_ok (const zoeDNA dna) {
	int OK = 1;
		
	if (dna->c5[0] + dna->c5[1] + dna->c5[2] + dna->c5[3] < 100) {
		zoeE("%s contains fewer than 100 unmasked nucleotides\n", dna->def);
		OK = 0;
	}
	if (dna->f5[4] > 0.95) {
		zoeE("%s contains more than 95%% masked nucleotides\n", dna->def);
		OK = 0;
	}

	return OK;
}

int gene_in_intron (const zoeCDS a, const zoeCDS b) {
	int        i;
	zoeFeature intron;
	
	/* assumed to overlap */
	for (i = 0; i < b->introns->size; i++) {
		intron = b->introns->elem[i];
		if (a->start > intron->start && a->end < intron->end) {
			return 1;
		}
	}
	
	return 0;
}

zoeFeatureVec get_xdef (const zoeDNA dna, const zoeFeatureTable ft, strand_t strand) {
	int           i;
	zoeFeature    f;
	zoeFeatureVec vec = zoeNewFeatureVec();
	
	/* gather features from correct strand */
	for (i = 0; i < ft->vec->size; i++) {
		f = ft->vec->elem[i];
		
		if (f->strand == strand) {
			zoePushFeatureVec(vec, f);
		} else if (f->strand == '=') {
			zoePushFeatureVec(vec, f);
			vec->last->strand = '+';
			zoePushFeatureVec(vec, f);
			vec->last->strand = '-';
		}
	}
	
	/* anti-features */
	for (i = 0; i < vec->size; i++) {
		f = vec->elem[i];
		if (f->strand == '-') zoeAntiFeature(f, dna->length);
	}
	
	return vec;
}

/* debugging */

static void debug_output (const zoeTrellis t) {
	int i, label;
	
	zoeO("TRELLIS\n");
	for (i = 0; i < t->dna->length; i++) {
		zoeO("%d", i);
		for (label = 0; label < zoeLABELS; label++) {
			if (t->score[label] == NULL) continue;
			if (t->score[label][i] == MIN_SCORE) zoeO("\t.:.:.");
			else zoeO("\t%d:%d", (int)t->score[label][i], t->trace[label][i]);
		}
		zoeO("\n");
	}
	
	zoeO("KEEP\n");
	for (i = 0; i < t->keep->size; i++) {
		zoeO("%d:%d\t", i, t->jump->elem[i]);
		zoeWriteFeature(stdout, t->keep->elem[i]);
	}
	
}


static void report_bug (zoeCDS gene) {
	int i;
	zoeE("unexpected gene errors found\n");
	zoeE("error list:\n");
	for (i = 0; i < gene->errors->size; i++) {
		zoeE("%s\n", gene->errors->elem[i]);
	}
	zoeE("warning list:\n");
	for (i = 0; i < gene->warnings->size; i++) {
		zoeE("%s\n", gene->warnings->elem[i]);
	}
	zoeE("exons:\n");
	for (i = 0; i < gene->exons->size; i++) {
		zoeWriteFeature(stderr, gene->exons->elem[i]);
	}
	
	zoeExit("please report bug including version, hmm, and sequence");
}

static void xdebug (const zoeTrellis t) {
	int i;
	char u0[16], u1[16], u2[16], a0[16], a1[16], a2[16];
	
	zoeO("xdebug strand %s\n", t->dna->def);
	for (i = 0; i < t->dna->length; i++) {
		if (t->scanner[Coding]->subscanner[0]->uscore) {
			zoeScore2Text(t->scanner[Coding]->subscanner[0]->uscore[i], u0);
			zoeScore2Text(t->scanner[Coding]->subscanner[1]->uscore[i], u1);
			zoeScore2Text(t->scanner[Coding]->subscanner[2]->uscore[i], u2);
		} else {
			strcpy(u0, ".");
			strcpy(u1, ".");
			strcpy(u2, ".");
		}
		if (t->scanner[Coding]->subscanner[0]->ascore) {
			zoeScore2Text(t->scanner[Coding]->subscanner[0]->ascore[i], a0);
			zoeScore2Text(t->scanner[Coding]->subscanner[1]->ascore[i], a1);
			zoeScore2Text(t->scanner[Coding]->subscanner[2]->ascore[i], a2);
		} else {
			strcpy(a0, ".");
			strcpy(a1, ".");
			strcpy(a2, ".");
		}
		zoeO("%d\t%s\t%s\t%s\t%s\t%s\t%s\n", i, u0, u1, u2, a0, a1, a2);
	}
	
	/*zoeExit("xdebug complete");*/
}

/* decoding */
		
zoeVec parse_strand (const zoeHMM hmm, const zoeDNA dna, const zoeFeatureTable ft, strand_t strand) {
	zoeTrellis    trellis;
	zoeVec        genes;
	zoeFeatureVec vec = NULL;
	
	if (zoeOption("-xdef")) vec = get_xdef(dna, ft, strand);
	trellis = zoeNewTrellis(dna, hmm, vec);
		
	genes = zoePredictGenes(trellis);
	if (zoeOption("-debug")) debug_output(trellis);
	if (zoeOption("-xdebug")) xdebug(trellis);
	
	zoeDeleteTrellis(trellis);
	if (vec) zoeDeleteFeatureVec(vec);
	
	return genes;
}

zoeVec parse_dna (const zoeHMM hmm, const zoeDNA plus_dna, zoeFeatureTable ft) {
	zoeDNA anti_dna;
	zoeVec plus_genes = NULL, anti_genes = NULL, genes, keep;
	zoeCDS gene, a, b;
	int    i, j, both_passed;
	char   id[64], name[256];
	
	
	/* plus */
	if (zoeOption("-plus")) {
		plus_genes = parse_strand(hmm, plus_dna, ft, '+');
		anti_genes = zoeNewVec();
	}
	else if (zoeOption("-minus")) {
		anti_dna = zoeAntiDNA(plus_dna->def, plus_dna);
		anti_genes = parse_strand(hmm, anti_dna, ft, '-');
		for (i = 0; i < anti_genes->size; i++) {
			zoeAntiCDS(anti_genes->elem[i], anti_dna->length);
		}
		zoeDeleteDNA(anti_dna);
		plus_genes = zoeNewVec();
	} else {
		plus_genes = parse_strand(hmm, plus_dna, ft, '+');
		anti_dna = zoeAntiDNA(plus_dna->def, plus_dna);
		anti_genes = parse_strand(hmm, anti_dna, ft, '-');
		for (i = 0; i < anti_genes->size; i++) {
			zoeAntiCDS(anti_genes->elem[i], anti_dna->length);
		}
		zoeDeleteDNA(anti_dna);
	}
	
	/****************************************\
		Both strands, sort out differences
	\****************************************/
	
	/* merge all genes into same vector, filter off tiny genes */
	genes = zoeNewVec();
	for (i = 0; i < plus_genes->size; i++) {
		gene = plus_genes->elem[i];
		
		if (!gene->OK) {
			if (strcmp(gene->errors->last, "gene:no_exons")) report_bug(gene);
			else continue;
		}
		
		if (gene->tx->length < SNAP_MIN_CDS) continue;
		if (gene->score < SNAP_MIN_SCORE) continue;
		
		zoePushVec(genes, gene);
	}
	for (i = 0; i < anti_genes->size; i++) {
		gene = anti_genes->elem[i];
		
		if (!gene->OK) {
			if (strcmp(gene->errors->last, "gene:no_exons")) report_bug(gene);
			else continue;
		}
		
		if (gene->tx->length < SNAP_MIN_CDS) continue;
		if (gene->score < SNAP_MIN_SCORE) continue;
		
		zoePushVec(genes, gene);
	}
	zoeDeleteVec(plus_genes); /* just containers, genes still kept */
	zoeDeleteVec(anti_genes);
	
	/* sort genes by score */
	qsort(genes->elem, genes->size, sizeof(zoeCDS), cmp_genes_ptr_by_score);

	/* compare genes */
	for (i = 0; i < genes->size; i++) {
		a = genes->elem[i];
		for (j = i+1; j < genes->size; j++) {
			b = genes->elem[j];
			if (!zoeCDSsOverlap(a, b)) continue;
			
			if (gene_in_intron(a, b) || gene_in_intron(b, a)) {
				if (zoeOption("-info")) zoeE("intronic genes %s %s\n", a->name, b->name);
				continue;
			}
			
			both_passed = 0;
			if      (a->score > SNAP_OVERLAP && b->score > SNAP_OVERLAP) both_passed = 1;
			else if (a->score > SNAP_OVERLAP)  b->score = MIN_SCORE;
			else if (b->score > SNAP_OVERLAP)  a->score = MIN_SCORE;
			else if (a->score > b->score) b->score = MIN_SCORE;
			else if (b->score > a->score) a->score = MIN_SCORE;
			else  {
				a->score = MIN_SCORE;
				b->score = MIN_SCORE;
			}
			
			/* may do something else when both pass */
			if (both_passed) {
				if (zoeOption("-info")) zoeE("overlap passed threshold %s %s\n", a->name, b->name);
			}
		}
	}
	
	/* throw out MIN_SCORE genes*/
	keep = zoeNewVec();
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		if (gene->score == MIN_SCORE) zoeDeleteCDS(gene);
		else                          zoePushVec(keep, gene);
	}
	zoeDeleteVec(genes);
	
	/* sort final genes by coordinate */
	qsort(keep->elem, keep->size, sizeof(zoeCDS), cmp_genes_ptr_by_start);
	
	/* rename genes */
	strncpy(id, plus_dna->def, 63);
	for (i = 0; i < strlen(plus_dna->def); i++) {
		if (isspace((int)plus_dna->def[i]) || i == 63) {
			id[i] = '\0';
			break;
		}
	}
	
	for (i = 0; i < keep->size; i++) {
		if (zoeOption("-name")) sprintf(name, "%s-%s.%d", id, zoeOption("-name"), i+1);
		else                    sprintf(name, "%s-snap.%d", id, i+1);
		
		gene = keep->elem[i];
		
		/* edit gene name */
		zoeFree(gene->name);
		gene->name = zoeMalloc(strlen(name) +1);
		strcpy(gene->name, name);
		
		/* edit transcript name */
		zoeFree(gene->tx->def);
		gene->tx->def = zoeMalloc(strlen(name) +1);
		strcpy(gene->tx->def, name);
		
		/* edit protein name */
		zoeFree(gene->aa->def);
		gene->aa->def = zoeMalloc(strlen(name) +1);
		strcpy(gene->aa->def, name);
		
		/* edit names in various feature vectors */
		edit_names(gene->exons,   name);
		edit_names(gene->introns, name);
		/*edit_names(gene->source,  name);*/
	}
	
	return keep;
}


/* output formatting */


void ace_output (const zoeDNA dna, const zoeVec genes) {
	int        i, j, start, end;
	zoeCDS     cds;
	zoeFeature exon;
	
	zoeO("\nSequence \"%s\"\n", dna->def);
	for (i = 0; i < genes->size; i++) {
		cds = genes->elem[i];
		
		if (cds->strand == '+') {
			start = cds->start +1;
			end   = cds->end   +1;
		} else {
			start = cds->end   +1;
			end   = cds->start +1;
		}
		
		zoeO("Subsequence %s %d %d\n", cds->name, start, end);
		
	}
	
	for (i = 0; i < genes->size; i++) {
		cds = genes->elem[i];
		zoeO("\nSequence %s\n", cds->name);
		if (cds->strand == '+') {
			for (j = 0; j < cds->exons->size; j++) {
				exon = cds->exons->elem[j];
				start = exon->start - cds->start +1;
				end   = exon->end   - cds->start +1;
				zoeO("Source_Exons %d %d\n", start, end);
			}
		} else {
			for (j = 0; j < cds->exons->size; j++) {
				exon = cds->exons->elem[j];
				end   = cds->end - exon->start +1;
				start = cds->end - exon->end   +1;
				zoeO("Source_Exons %d %d\n", start, end);
			}
		}
			zoeO("Method SNAP\n");
	}
}

void gff_output (const zoeDNA dna, const zoeVec genes) {
	int i, j;
	zoeCDS gene;

	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		for (j = 0; j < gene->exons->size; j++) {
			zoeWriteGFF(stdout, gene->exons->elem[j], "SNAP", dna->def);
		}
	}
}

void zoe_output (const zoeDNA dna, const zoeVec genes) {
	int           i, j;
	zoeCDS        gene;
	zoeFeatureVec fv;
		
	zoeO(">%s\n", dna->def);
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		if (zoeOption("-extra")) fv = gene->source;
		else                     fv = gene->exons;
		for (j = 0; j < fv->size; j++) {
			zoeWriteFeature(stdout, fv->elem[j]);
		}	
	}
}

void help (void) {

zoeM(stdout, 35,

"The general form of the snap command line is:\n",

"    snap <HMM file> <FASTA file> [options]\n",

"HMM file:\n",

"    The most convenient way to specify the HMM file is by name. This requires",
"    that the ZOE environment variable is set. In this case, snap will look",
"    for the HMM file in $ZOE/HMM. You may also specify the HMM file by an",
"    explicit path. The following are equivalent if $ZOE is in /usr/local:\n",
"        snap C.elegans.hmm ...",
"        snap /usr/local/Zoe/HMM/C.elegans.hmm ...",
"        snap worm ...  # there are a few convenient aliases in $ZOE/HMM\n",

"FASTA file:\n",

"    If you have several sequences to analyze, it is more efficient to run",
"    snap on a concatenated FASTA file rather than separate runs on single",
"    sequence files. The seqeuence may be in a compressed format\n",

"    If sequences have been masked with lowercase letters, use -lcmask to",
"    prevent exons from appearing in masked DNA.\n",

"Output:\n",

"    Annotation is reported to stdout in a non-standard format (ZFF). You can",
"    change to GFF or ACEDB with the -gff or -ace options. Proteins and",
"    transcripts are reported to FASTA files with the -aa and -tx options.\n",

"External definitions:\n",

"    SNAP allows you to adjust the score of any sequence model at any point",
"    in a sequence. This behavior is invoked by giving a ZFF file to SNAP:\n",

"        snap <hmm> <sequence> -xdef <ZFF file>\n",

"    Each feature description uses the 'group' field to issue a command:\n",

"        SET     set the score",
"        ADJ     adjust the score up or down",
"        OK      set non-cannonical scores\n",

"     >FOO",
"     Acceptor 120 120 + +50 . . . SET  (sets an Acceptor to 50)",
"     Donor    212 212 + -20 . . . ADJ  (lowers a Donor by -20)",
"     Inter    338 579 +  -2 . . . ADJ  (lowers Inter by -2 in a range)",
"     Coding   440 512 -  +3 . . . ADJ  (raises Coding by +3 in a range)",
"     Donor    625 638 +  -5 . . . OK   (sets range of odd Donors to -5)\n",

"If the output has scrolled off your screen, try 'snap -help | more'"

);

exit(0);

}

