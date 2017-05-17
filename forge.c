/*****************************************************************************\
 forge.c

Parameter estimation program for SNAP

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

/* model struct */
struct model {
	zoeLabel   label;
	int        length;
	int        order;
	zoeModel   counts;
	zoeModel   scores;
};

/* prototypes */
void help (void);
float gc_content (void);
float repeat_content (void);
void openData (void);
void closeData (void);
int getData (void);
zoeFeatureVec extractAcceptor (void);
zoeFeatureVec extractDonor (void);
zoeFeatureVec extractStart (void);
zoeFeatureVec extractStop (void);
zoeFeatureVec extractInter (void);
zoeFeatureVec extractIntron (void);
zoeFeatureVec extractCoding (void);
zoeFeatureVec extractUTR3 (void);
zoeFeatureVec extractUTR5 (void);
zoeFeatureVec extractPolyA (void);

void init_boosting (void);
zoeVec define_models (void);
void countModels (void);
void createModel (zoeModel, const zoeModel);
void processWMM (zoeModel, const zoeModel);
void processLUT (zoeModel, const zoeModel);
void createModels (void);
void outputModels (void);
void outputPolyA (void);

void countDurations (void);
void outputDurations (void);

void countTransitions (void);
void outputTransitions (void);

void calcPhasePref (void);
void outputPhasePref (void);

void buildHMMs (void);

/* global: DNA & Features */
char * DNA_file;
char * ANN_file;
zoeFile DNA_stream;
zoeFile ANN_stream;
zoeDNA DNA;           /* the current DNA */
zoeDNA ANTI;          /* reverse complement of DNA */
zoeFeatureTable ANN;  /* the current annotation */
zoeCDS GENE;          /* the current GENE */
zoeFeatureVec Features[zoeLABELS];
int UTR5_Length = 50;
int UTR3_Length = 50;
int UTR5_Offset = 10;
int UTR3_Offset = 10;
int POLYA_LOOK  = 500;

/* global: phases, transitions, durations */
int PhasePref[18];
int Transition[zoeLABELS][zoeLABELS];
#define DURATION_LIMIT 5000
float IntronD[DURATION_LIMIT];
float InterD[DURATION_LIMIT];
float ExonD[DURATION_LIMIT];
float EinitD[DURATION_LIMIT];
float EtermD[DURATION_LIMIT];
float EsnglD[DURATION_LIMIT];

/* global: models and scores */
zoeVec Models;

/* global: Esngl static length data */
static int EsnglStaticData[96] = {
	201, 216, 222, 234, 246, 261, 273, 279, 291, 294, 
	306, 312, 312, 324, 330, 339, 345, 354, 363, 372, 
	384, 384, 384, 387, 399, 408, 411, 420, 432, 441, 
	453, 468, 480, 489, 504, 516, 528, 540, 555, 567, 
	582, 594, 609, 618, 633, 645, 666, 681, 702, 720, 
	735, 747, 762, 777, 792, 813, 831, 849, 873, 891, 
	909, 933, 945, 969, 987, 1011, 1020, 1020, 1032, 1044, 
	1080, 1095, 1140, 1173, 1200, 1236, 1272, 1326, 1350, 1401, 
	1428, 1449, 1506, 1557, 1602, 1662, 1773, 1830, 1914, 2031, 
	2112, 2226, 2418, 2637, 2913
};


/* global: min counts */
static int MIN_COUNTS = 0;
static int MIN_INTRON = 30;
static int VERBOSE = 0;
static zoeHash WEIGHT = NULL; /* boosting */

static char usage[]  = "\n\
FORGE - training program for SNAP (version 2006-07-28)\n\n\
usage: forge [options] <ann> <dna> [options]\n\
options:\n\
  -help\n\
  -verbose\n\
  -pseudocount <float>  [1]   (absolute number for all models)\n\
  -pseudoCoding <float> [0.0] (eg. 0.05)\n\
  -pseudoIntron <float> [0.0]\n\
  -pseudoInter <float>  [0.0]\n\
  -min-counts           [0]\n\
  -lcmask [-fragmentN]\n\
  -utr5-length <int>    [50]\n\
  -utr5-offset <int>    [10]\n\
  -utr3-length <int>    [50]\n\
  -utr3-offset <int>    [10]\n\
  -explicit <int>       [250]\n\
  -min-intron <int>     [30]\n\
  -boost <file>  (file of ID <int>)\n\
";

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char **argv) {
	int i, label;
	
	/* preamble */
	zoeSetProgramName(argv[0]);
	zoeSetOption("-help", 0);
	zoeSetOption("-verbose", 0);
	zoeSetOption("-pseudocount", 1);
	zoeSetOption("-pseudoInter", 1);
	zoeSetOption("-pseudoIntron", 1);
	zoeSetOption("-pseudoCoding", 1);
	zoeSetOption("-min-counts", 1);
	zoeSetOption("-lcmask", 0);
	zoeSetOption("-fragmentN", 0);
	zoeSetOption("-inter+intron", 0);
	zoeSetOption("-polyA", 1);
	zoeSetOption("-utr5-length", 1);
	zoeSetOption("-utr3-length", 1);
	zoeSetOption("-utr5-offset", 1);
	zoeSetOption("-utr3-offset", 1);
	zoeSetOption("-explicit", 1);
	zoeSetOption("-min-intron", 1);
	zoeSetOption("-boost", 1);
	zoeParseOptions(&argc, argv);
	
	if (zoeOption("-help")) help();
	
	if (argc != 3) {
		zoeE("%s", usage);
		exit(1);
	}
	
	ANN_file = argv[1];
	DNA_file = argv[2];
	
	/* Minimum counts */
	if (zoeOption("-min-counts")) MIN_COUNTS = atoi(zoeOption("-min-counts"));
	if (zoeOption("-verbose")) VERBOSE = 1;
	
	/* UTR lengths and PolyA scan */
	if (zoeOption("-utr5-length")) UTR5_Length = atoi(zoeOption("-utr5-length"));
	if (zoeOption("-utr5-offset")) UTR5_Offset = atoi(zoeOption("-utr5-offset"));
	if (zoeOption("-utr3-length")) UTR3_Length = atoi(zoeOption("-utr3-length"));
	if (zoeOption("-utr3-offset")) UTR3_Offset = atoi(zoeOption("-utr3-offset"));
	if (zoeOption("-polyA")) POLYA_LOOK = atoi(zoeOption("-polyA"));
	
	/* initialze global data structures */
	DNA        = NULL;
	ANTI       = NULL;
	ANN        = NULL;
	GENE       = NULL;
	for (label = 0; label < zoeLABELS; label++) {
		Features[label] = NULL;
		for (i = 0; i < zoeLABELS; i++) Transition[label][i] = 0;
	}
	for (i = 0; i < DURATION_LIMIT; i++) {
		IntronD[i] = 0; InterD[i]  = 0; ExonD[i]  = 0;
		EinitD[i]  = 0; EtermD[i]  = 0; EsnglD[i] = 0;
	}
	for (i = 0; i < 18; i++) PhasePref[i] = 0;
	Models = define_models();
	
	/* Others */
	if (zoeOption("-min-intron")) MIN_INTRON = atoi(zoeOption("-min-intron"));
	if (zoeOption("-boost")) init_boosting();
	
	/* data loop */
	if (VERBOSE) zoeE("counting");
	openData();
	while (getData()) {
		calcPhasePref();
		countModels();
		countDurations();
		countTransitions();
		if (VERBOSE) zoeE(".");
	}
	closeData();
	if (VERBOSE) zoeE("done\n");
	
	/* models */
	createModels();
	
	/* output components */	
	outputPhasePref();
	outputDurations();
	outputTransitions();
	outputModels();
	
	return 0;
}

/*****************************************************************************\
 model functions
\*****************************************************************************/

/*
void calcMeanModelScore (void) {
	float total, var;
	int i, j;
	struct model * m;
	
	for (i = 0; i < Models->size; i++) {
		m = Models->elem[i];
		
		total = 0;
		for (j = 0; j < m->scores->size; j++) total += m->scores->elem[j] - SCORE_SHELF;
		m->mean_score = zoeDivide(total, m->scores->size);

		total = 0;
		for (j = 0; j < m->scores->size; j++)
			total += pow(m->scores->elem[j] - m->mean_score - SCORE_SHELF, 2);
		var = zoeDivide(total, (m->scores->size - 1));
		m->std_dev = sqrt(var);		
	}
}
*/

void init_boosting (void) {
	char    id[64];
	int     weight;
	size_t  set;
	zoeFile file;
	
	WEIGHT = zoeNewHash();
	
	file = zoeOpenFile(zoeOption("-boost"));
	while (fscanf(file.stream, "%s %d", id, &weight) != EOF) {
		set = weight;
		zoeSetHash(WEIGHT, id, (void*)set);
	}
	
	zoeCloseFile(file);
}

void countModels (void) {
	int i, j, k, w;
	char * id;
	struct model * m;
	zoeCounter c;
	zoeFeature f;
	
	for (i = 0; i < Models->size; i++) {
		m = Models->elem[i];
		c = zoeNewCounter(DNA, ANTI, m->counts);
		switch (m->label) {
			case Acceptor: case Donor: case Start: case Stop: case PolyA:
				for (j = 0; j < Features[m->label]->size; j++) {
					f = Features[m->label]->elem[j];
					if (WEIGHT) {
						id = DNA->def;
						w = (size_t)zoeGetHash(WEIGHT, id);
						for (k = 0; k < w; k++) c->count(c, f->start);
					} else c->count(c, f->start);
				}
				break;
			case Inter: case Intron: case Coding: case UTR5: case UTR3:
				for (j = 0; j < Features[m->label]->size; j++) {
					f = Features[m->label]->elem[j];
					if (WEIGHT) {
						id = DNA->def;
						w = (size_t)zoeGetHash(WEIGHT, id);
						for (k = 0; k < w; k++) c->countf(c, f);
					} else c->countf(c, f);
				}
				break;
			default: zoeExit("countModels not handled yet");
		}
		zoeDeleteCounter(c);
	}
}

void outputModels (void) {
	int i;
	struct model * m;
	FILE * file;
	char filename[256];
	
	if (VERBOSE) zoeE("writing models");
	for (i = 0; i < Models->size; i++) {
		if (VERBOSE) zoeE(".");
		m = Models->elem[i];
		
		/* counts */
		sprintf(filename, "%s-%d-%d.count", m->counts->name, m->order, m->length);
		file = fopen(filename, "w");
		zoeWriteModel(file, m->counts);
		fclose(file);
		
		/* scores */
		sprintf(filename, "%s-%d-%d.model", m->counts->name, m->order, m->length);
		file = fopen(filename, "w");
		zoeWriteModel(file, m->scores);
		fclose(file);
		
		
	}
	if (VERBOSE) zoeE("done\n");
}

static void fractional_pseudocounts (zoeModel counts, float fraction) {
	int i, slots, total = 0;
	float pseudo;
	
	slots = zoePOWER[counts->symbols][counts->length];
	for (i = 0; i < slots; i++) total += counts->data[i];
	pseudo = fraction * total / slots;
	for (i = 0; i < slots; i++) counts->data[i] += pseudo;
}

void createModels (void) {
	int i;
	struct model * m;
	
	if (VERBOSE) zoeE("creating models");
	for (i = 0; i < Models->size; i++) {
		if (VERBOSE) zoeE(".");
		m = Models->elem[i];
		
		zoeDeambiguateModel(m->counts);
		
		if (m->label == Inter  && zoeOption("-pseudoInter")) {
			fractional_pseudocounts(m->counts, atof(zoeOption("-pseudoInter")));
		} else if (m->label == Intron && zoeOption("-pseudoIntron")) {
			fractional_pseudocounts(m->counts, atof(zoeOption("-pseudoIntron")));
		} else if (m->label == Coding && zoeOption("-pseudoCoding")) {
			fractional_pseudocounts(m->counts->submodel[0], atof(zoeOption("-pseudoCoding")));
			fractional_pseudocounts(m->counts->submodel[1], atof(zoeOption("-pseudoCoding")));
			fractional_pseudocounts(m->counts->submodel[2], atof(zoeOption("-pseudoCoding")));

		}
		createModel(m->scores, m->counts);
	}
	if (VERBOSE) zoeE("done\n");
}

void createModel (zoeModel scores, const zoeModel counts) {
	int i;
	
	switch (scores->type) {
		case SAM: case SDT: case CDS:
			for (i = 0; i < scores->submodels; i++)
				createModel(scores->submodel[i], counts->submodel[i]);
			break;
		case WMM: processWMM(scores, counts); break;
		case LUT: processLUT(scores, counts); break;
		case TRM: break;
		default:  zoeExit("createModel unable to handle %s\n", scores->name);
	}
}

void processWMM (zoeModel scores, const zoeModel counts) {
	int     i, j;
	float   total, frac;
	score_t score;
	int     index;
	
	/* min counts */
	for (i = 0; i < 4; i++) {
		if (counts->data[i] < MIN_COUNTS) counts->data[i] = 0;
	}
			
	/* compute total counts */
	total = 0;
	for (i = 0; i < 4; i++) total += counts->data[i];
	
	/* create scores */
	for (i = 0; i < scores->length; i++) {
		for (j = 0; j < 4; j++) {
			index  = i * 4 + j;
			if (total <= 0 || counts->data[index] == 0) {
				score = MIN_SCORE;
			} else {
				frac  = counts->data[index] / total;
				score = zoeLog2(frac/0.25);
			}
			scores->data[index] = score;
		}
	}
}

void processLUT (zoeModel scores, const zoeModel counts) {
	int     i, j;
	int     p = zoePOWER[4][counts->length];
	float   total, frac;
	score_t score;
		
	for (i = 0; i < p; i+= 4) {
		
		/* min counts */
		for (j = 0; j < 4; j++) {
			if (counts->data[i+j] < MIN_COUNTS) counts->data[i+j] = 0;
		}
	
		/* get total counts for this context */
		total = 0;
		for (j = 0; j < 4; j++) total += counts->data[i+j];
				
		/* create scores for this context */
		for (j = 0; j < 4; j++ ) {
			if (total <= 0 || counts->data[i+j] == 0) {
				score = MIN_SCORE;
			} else {
				frac  = counts->data[i+j] / total;
				score = zoeLog2(frac/0.25);
			}
			scores->data[i+j] = score;
		}
	}
}


/*****************************************************************************\
 transition functions
\*****************************************************************************/

void countTransitions (void) {
	int i;
	
	for (i = 0; i < Features[Coding]->size; i++) {
		if (Features[Coding]->elem[i]->label == 0) continue;
		switch (Features[Coding]->elem[i]->label) {
			case Einit: Transition[Inter][Einit]++;  break;
			case Esngl: Transition[Inter][Esngl]++;  break;
			case Eterm: Transition[Intron][Eterm]++; break;
			case Exon:  Transition[Intron][Exon]++;  break;
			default: break;
		}
	}
}

void outputTransitions (void) {
	FILE * file;
	float total;
	
	file = fopen("transitions", "w");
	
	total = Transition[Inter][Einit] + Transition[Inter][Esngl];
	zoeS(file, "Inter Einit %f\n", Transition[Inter][Einit]/total);
	zoeS(file, "Inter Esngl %f\n", Transition[Inter][Esngl]/total);
	
	total = Transition[Intron][Exon] + Transition[Intron][Eterm];
	zoeS(file, "Intron Exon %f\n", Transition[Intron][Exon]/total);
	zoeS(file, "Intron Eterm %f\n", Transition[Intron][Eterm]/total);
	
	fclose(file);
}


/*****************************************************************************\
 duration functions
\*****************************************************************************/

void countDurations (void) {
	int i, length;
	zoeFeature f;
	
	for (i = 0; i < Features[Coding]->size; i++) {
		f = Features[Coding]->elem[i];
		if (f->label == 0) continue;
		length = f->end - f->start +1;
		if (length >= DURATION_LIMIT) length = DURATION_LIMIT -1;
		switch (f->label) {
			case Exon:  ExonD[length]++;  break;
			case Einit: EinitD[length]++; break;
			case Eterm: EtermD[length]++; break;
			case Esngl: EsnglD[length]++; break;
			default: zoeExit("pollution of exons");
		}
	}
	
	for (i = 0; i < Features[Intron]->size; i++) {
		f = Features[Intron]->elem[i];
		if (f->label == 0) continue;
		length = f->end - f->start +1;
		if (length >= DURATION_LIMIT) length = DURATION_LIMIT -1;
		IntronD[length]++;
	}
	
	for (i = 0; i < Features[Inter]->size; i++) {
		f = Features[Inter]->elem[i];
		if (f->label == 0) continue;
		length = f->end - f->start +1;
		if (length >= DURATION_LIMIT) length = DURATION_LIMIT -1;
		InterD[length]++;
	}
	
}

static void outputDuration (zoeLabel label, int limit) {
	char  * name;
	float * ary;
	FILE  * file;
	char    filename[1024];
	char    s[64];
	float   total, length, mean;
	int     i, j;
	
	switch (label) {
		case Einit:  name = "Einit";  ary = EinitD;  break;
		case Eterm:  name = "Eterm";  ary = EtermD;  break;
		case Esngl:  name = "Esngl";  ary = EsnglD;  break;
		case Exon:   name = "Exon";   ary = ExonD;   break;
		case Intron: name = "Intron"; ary = IntronD; break;
		case Inter:  name = "Inter";  ary = InterD;  break;
		default:
			name = ""; /* shush */
			ary = NULL;
			zoeExit("outputDuration does not support %d", label);
	}
	
	/* explicit duration */
	sprintf(filename, "%s-explicit.duration", name);
	file = fopen(filename, "w");
	
	zoeS(file, "%s 2\n", name);
	zoeS(file, "\tDEFINED 0 %d\n", limit -1);
	
	total = 0;
	length = 0;
	for (i = 0; i < DURATION_LIMIT; i++) {
		total += ary[i];
		length += ary[i] * i;
	}
	mean = length / total;
	
	for (i = 0; i < limit; i+=5) {
		zoeS(file, "\t\t");
		for (j = 0; j < 5; j++) {
			if (label == Intron && (i+j) < MIN_INTRON) {
				zoeS(file, ". ");
			} else {
				zoeScore2Text(zoeFloat2Score(ary[i+j]/total), s);
				zoeS(file, "%s ", s);
			}
		}
		zoeS(file, "\n");
	}
	
	zoeS(file, "\tGEOMETRIC %d -1\n", limit);
	zoeS(file, "\t\t%d\n", (int)mean);
	
	fclose(file);
	
	/* geometric duration */
	sprintf(filename, "%s-geometric.duration", name);
	file = fopen(filename, "w");
	
	zoeS(file, "%s 1\n", name);
	zoeS(file, "\tGEOMETRIC 0 -1\n");
	zoeS(file, "\t\t%d\n", (int)mean);
	fclose(file);

}

static float ave_counts (float * count, int pos, int window) {
	int i;
	float t, mean;
	
	t = 0;
	for (i = pos - window; i <= pos + window; i++) t += count[i];
	mean = t / (float)(window * 2 + 1);
	return mean;
}

static void smoothCounts (float * count) {
	float   tmp[DURATION_LIMIT];
	int     i, j;
	float   a, m;
	zoeFVec s1;
	int   window, step;
	
	for (i = 0; i < DURATION_LIMIT; i++) tmp[i] = 0;
	s1 = zoeNewFVec();
	
	window = 20;
	step   = 10;
	for (i = step; i < DURATION_LIMIT - window; i += step) {
		a = ave_counts(count, i, window);
		zoePushFVec(s1, a);
	}
	
	for (i = 0; i < s1->size -1; i++) {
		m = (s1->elem[i+1] - s1->elem[i]) / step;
		for (j = 0; j < step; j++) {
			tmp[i*step +j] = s1->elem[i] + m * j;
		}
	}
	
	zoeDeleteFVec(s1);
	
	/* overwrite count array */
	for (i = 0; i < DURATION_LIMIT; i++) count[i] = tmp[i] + 0.0001; /* prevent 0 errors */
}

static void addFixedEsnglCounts (void) {
	int i;
	for (i = 0; i < 96; i++) EsnglD[EsnglStaticData[i]]++;
}

void outputDurations (void) {
	int fixed = 250;
	
	if (zoeOption("-explicit")) fixed = atoi(zoeOption("-explicit"));
	
	if (VERBOSE) zoeE("durations...");
	
	smoothCounts(EinitD);
	outputDuration(Einit, fixed);
	
	smoothCounts(ExonD);
	outputDuration(Exon, fixed);
	
	smoothCounts(EtermD);
	outputDuration(Eterm, fixed);
	
	addFixedEsnglCounts();
	smoothCounts(EsnglD);
	outputDuration(Esngl, 3000);
	
	smoothCounts(IntronD);
	outputDuration(Intron, fixed);

	smoothCounts(InterD);
	outputDuration(Inter, fixed);
	
}


/*****************************************************************************\
 phasepref functions
\*****************************************************************************/

void calcPhasePref (void) {
	int i;
	zoeFeature prev, exon;
	int index = -1; /* impossible value */
		
	for (i = 1; i < GENE->exons->size; i++) {
		prev = GENE->exons->elem[i-1];
		exon = GENE->exons->elem[i];
	
		/* determine proper index for exon->intron transition */
		if (prev->label == Einit) {
			switch (prev->inc3) {
				case 0: index = 0; break; /* initial exon to phase 0 intron */
				case 1: index = 1; break; /* initial exon to phase 1 intron */
				case 2: index = 2; break; /* initial exon to phase 2 intron */
				default: zoeExit("impossible error1");
			}
		} else if (prev->label == Exon) {
			if (prev->inc3 == 0) {
				switch (exon->inc3) {
					case 0: index = 3; break; /* 0....0 exon */
					case 1: index = 4; break; /* 0....1 exon */
					case 2: index = 5; break; /* 0....2 exon */
					default: zoeExit("impossible error2");
				}
			} else if (prev->inc3 == 1) {
				switch (exon->inc3) {
					case 0: index = 6; break; /* 1....0 exon */
					case 1: index = 7; break; /* 1....1 exon */
					case 2: index = 8; break; /* 1....2 exon */
					default: zoeExit("impossible error3");
				}
			} else if (prev->inc3 == 2) {
				switch (exon->inc3) {
					case 0: index = 9;  break; /* 2....0 exon */
					case 1: index = 10; break; /* 2....1 exon */
					case 2: index = 11; break; /* 2....2 exon */
					default: zoeExit("impossible error4");
				}
			}
		}
		PhasePref[index]++;
	
		/* determine proper index for intron->exon transition */
		if (prev->inc3 == 0) {
			switch (exon->label) {
				case Exon:  index = 12; break; /* phase 0 intron to Exon */
				case Eterm: index = 13; break; /* phase 0 intron to Eterm */
				default: zoeExit("impossible error5 %s", GENE->name);
			}
		} else if (prev->inc3 == 1) {
			switch (exon->label) {
				case Exon:  index = 14; break; /* phase 1 intron to Exon */
				case Eterm: index = 15; break; /* phase 1 intron to Eterm */
				default: zoeExit("impossible error6");
			}
		} else if (prev->inc3 == 2) {
			switch (exon->label) {
				case Exon:  index = 16; break; /* phase 2 intron to Exon */
				case Eterm: index = 17; break; /* phase 2 intron to Eterm */
				default: zoeExit("impossible error7");
			}
		}
		PhasePref[index]++;
	}
}


void outputPhasePref (void) {
	int total;
	FILE * outfile;
	
	outfile = fopen("phaseprefs", "w");
	
	/* Einit to intron */
	total = PhasePref[0] + PhasePref[1] + PhasePref[2];
	zoeS(outfile, "%f\n%f\n%f\n",
		zoeDivide((float)PhasePref[0],(float)total),
		zoeDivide((float)PhasePref[1],(float)total),
		zoeDivide((float)PhasePref[2],(float)total));
	
	/* phase 0 Exon to intron */
	total = PhasePref[3] + PhasePref[4] + PhasePref[5];
	zoeS(outfile, "%f\n%f\n%f\n",
		zoeDivide((float)PhasePref[3],(float)total),
		zoeDivide((float)PhasePref[4],(float)total),
		zoeDivide((float)PhasePref[5],(float)total));

	/* phase 1 Exon to intron */
	total = PhasePref[6] + PhasePref[7] + PhasePref[8];
	zoeS(outfile, "%f\n%f\n%f\n",
		zoeDivide((float)PhasePref[6],(float)total),
		zoeDivide((float)PhasePref[7],(float)total),
		zoeDivide((float)PhasePref[8],(float)total));

	/* phase 2 Exon to intron */
	total = PhasePref[9] + PhasePref[10] + PhasePref[11];
	zoeS(outfile, "%f\n%f\n%f\n",
		zoeDivide((float)PhasePref[9],(float)total),
		zoeDivide((float)PhasePref[10],(float)total),
		zoeDivide((float)PhasePref[11],(float)total));

	/* phase 0 intron to Exon and Eterm */
	total = PhasePref[12] + PhasePref[13];
	zoeS(outfile, "%f\n%f\n",
		zoeDivide((float)PhasePref[12],(float)total),
		zoeDivide((float)PhasePref[13],(float)total));
	
	/* phase 1 intron to Exon and Eterm */
	total = PhasePref[14] + PhasePref[15];
	zoeS(outfile, "%f\n%f\n",
		zoeDivide((float)PhasePref[14],(float)total),
		zoeDivide((float)PhasePref[15],(float)total));
	
	/* phase 2 intron to Exon and Eterm */
	total = PhasePref[16] + PhasePref[17];
	zoeS(outfile, "%f\n%f\n",
		zoeDivide((float)PhasePref[16],(float)total),
		zoeDivide((float)PhasePref[17],(float)total));
	
	fclose(outfile);
}

/*****************************************************************************\
 extract
\*****************************************************************************/

zoeFeatureVec extractInter (void) {
	zoeFeature    inter;
	zoeFeatureVec vec = zoeNewFeatureVec();
	int           i;
	
	/* 5' intergenic */
	if (GENE->start > 0) {
		inter = zoeNewTriteFeature(Inter, 0, GENE->start -1, NULL);
		zoePushFeatureVec(vec, inter);
		zoeDeleteFeature(inter);
	}
	
	/* 3' intergenic */
	if (GENE->end < DNA->length -1) {
		inter = zoeNewTriteFeature(Inter, GENE->end +1, DNA->length -1, NULL);
 		zoePushFeatureVec(vec, inter);
		zoeDeleteFeature(inter);
	}
	
	/* merge intron and intergenic */
	if (zoeOption("-inter+intron")) {
		for (i = 0; i < GENE->introns->size; i++) {
			zoePushFeatureVec(vec, GENE->introns->elem[i]);
		}
	}
	return vec;
}

zoeFeatureVec extractAcceptor (void) {
	int           i, coor;
	zoeFeature    exon;
	zoeFeature    acc = zoeNewTriteFeature(Acceptor, 0, 1, NULL);
	zoeFeatureVec vec = zoeNewFeatureVec();
		
	for (i = 0; i < GENE->exons->size; i++) {
		exon = GENE->exons->elem[i];
		if ( !(exon->label == Exon || exon->label == Eterm) ) continue;
		coor = exon->start -1;
		acc->start = coor;
		acc->end   = coor;
		zoePushFeatureVec(vec, acc);
	}
	zoeDeleteFeature(acc);

	return vec;
}

zoeFeatureVec extractDonor (void) {
	int           i, coor;
	zoeFeature    exon;
	zoeFeatureVec vec = zoeNewFeatureVec();
	zoeFeature    don = zoeNewTriteFeature(Donor, 0, 0, NULL);
		
	for (i = 0; i < GENE->exons->size; i++) {
		exon = GENE->exons->elem[i];
		if ( !(exon->label == Exon || exon->label == Einit) ) continue;
		coor = exon->end +1;
		don->start = coor;
		don->end   = coor;
		zoePushFeatureVec(vec, don);
	}
	zoeDeleteFeature(don);

	return vec;
}

zoeFeatureVec extractStart (void) {
	int           i, coor;
	zoeFeature    start, exon;
	zoeFeatureVec vec = zoeNewFeatureVec();
		
	for (i = 0; i < GENE->exons->size; i++) {
		exon = GENE->exons->elem[i];
		if ( !(exon->label == Einit || exon->label == Esngl) ) continue;
		coor = exon->start;
		start = zoeNewTriteFeature(Start, coor, coor, NULL);
		zoePushFeatureVec(vec, start);
		zoeDeleteFeature(start);
	}

	return vec;
}

zoeFeatureVec extractStop (void) {
	int           i, coor;
	zoeFeature    stop, exon;
	zoeFeatureVec vec = zoeNewFeatureVec();
		
	for (i = 0; i < GENE->exons->size; i++) {
		exon = GENE->exons->elem[i];
		if ( !(exon->label == Eterm || exon->label == Esngl) ) continue;
		coor = exon->end -2;
		stop = zoeNewTriteFeature(Stop, coor, coor, NULL);
		zoePushFeatureVec(vec, stop);
		zoeDeleteFeature(stop);
	}
	
	return vec;
}

zoeFeatureVec extractIntron (void) {
	int           i;
	zoeFeatureVec vec = zoeNewFeatureVec();
		
	for (i = 0; i < GENE->introns->size; i++) zoePushFeatureVec(vec, GENE->introns->elem[i]);

	return vec;
}

zoeFeatureVec extractCoding (void) {
	int           i;
	zoeFeatureVec vec = zoeNewFeatureVec();
		
	for (i = 0; i < GENE->exons->size; i++) zoePushFeatureVec(vec, GENE->exons->elem[i]);
	
	return vec;
}

zoeFeatureVec extractUTR5 (void) {
	int           a, b;
	zoeFeature    utr;
	zoeFeatureVec vec = zoeNewFeatureVec();
		
	a = GENE->start - UTR5_Length - UTR5_Offset;
	b = GENE->start - UTR5_Length;
	
	if (a < 0) a = 0;
	if (b < 0) b = 0;
	
	utr = zoeNewTriteFeature(UTR5, a, b, NULL);
	zoePushFeatureVec(vec, utr);
	zoeDeleteFeature(utr);
	
	return vec;
}

zoeFeatureVec extractUTR3 (void) {
	int           a, b;
	zoeFeature    utr;
	zoeFeatureVec vec = zoeNewFeatureVec();
	
	a = GENE->end + UTR3_Offset;
	b = GENE->end + UTR3_Offset + UTR3_Length;
	
	if (a >= DNA->length) a = DNA->length -1;
	if (b >= DNA->length) b = DNA->length -1;
	
	utr = zoeNewTriteFeature(UTR3, a, b, NULL);
	zoePushFeatureVec(vec, utr);
	zoeDeleteFeature(utr);
	
	return vec;
}

zoeFeatureVec extractPolyA (void) {
	int           i, a, b;
	zoeFeature    polyA;
	zoeFeatureVec vec = zoeNewFeatureVec();
	score_t       score, max_score;
	coor_t        max_pos;
		
	/*
	
	Looking for a reasonable match to AATAAA downstream of stop codon.
	Data for hard-coded scoring matrix adapted from Tom Blumenthal.
	
	This is a hack and not a rigorous way to annotate the polyA signal.
	
	*/
	
	a = GENE->end + UTR3_Offset;
	b = GENE->end + POLYA_LOOK;
	if (b >= DNA->length) b = DNA->length;
	max_score = 0;
	max_pos = 0;
	
	for (i = a; i < b -5; i++) {
		score = 1;
	
		if      (DNA->s5[i] == 0) score *= 0.80;
		else if (DNA->s5[i] == 1) score *= 0.05;
		else if (DNA->s5[i] == 2) score *= 0.05;
		else if (DNA->s5[i] == 3) score *= 0.10;
		else if (DNA->s5[i] == 4) score *= 0;
	
		if      (DNA->s5[i+1] == 0) score *= 0.95;
		else if (DNA->s5[i+1] == 1) score *= 0.01;
		else if (DNA->s5[i+1] == 2) score *= 0.02;
		else if (DNA->s5[i+1] == 3) score *= 0.12;
		else if (DNA->s5[i+1] == 4) score *= 0;
	
		if      (DNA->s5[i+2] == 0) score *= 0.03;
		else if (DNA->s5[i+2] == 1) score *= 0.01;
		else if (DNA->s5[i+2] == 2) score *= 0.01;
		else if (DNA->s5[i+2] == 3) score *= 0.95;
		else if (DNA->s5[i+2] == 4) score *= 0;
	
		if      (DNA->s5[i+3] == 0) score *= 0.85;
		else if (DNA->s5[i+3] == 1) score *= 0.01;
		else if (DNA->s5[i+3] == 2) score *= 0.13;
		else if (DNA->s5[i+3] == 3) score *= 0.01;
		else if (DNA->s5[i+3] == 4) score *= 0;
	
		if      (DNA->s5[i+4] == 0) score *= 0.96;
		else if (DNA->s5[i+4] == 1) score *= 0.02;
		else if (DNA->s5[i+4] == 2) score *= 0.01;
		else if (DNA->s5[i+4] == 3) score *= 0.01;
		else if (DNA->s5[i+4] == 4) score *= 0;
	
		if      (DNA->s5[i+5] == 0) score *= 0.96;
		else if (DNA->s5[i+5] == 1) score *= 0.01;
		else if (DNA->s5[i+5] == 2) score *= 0.01;
		else if (DNA->s5[i+5] == 3) score *= 0.02;
		else if (DNA->s5[i+5] == 4) score *= 0;
	
		if (score > max_score) {
			max_score = score;
			max_pos   = i;
		}
	}
	
	if (max_score > 0.01) {
		polyA = zoeNewTriteFeature(PolyA, max_pos, max_pos, NULL);
		zoePushFeatureVec(vec, polyA);
		zoeDeleteFeature(polyA);
	}

	return vec;
}


/*****************************************************************************\
 misc functions
\*****************************************************************************/

void help (void) {
	zoeO("sorry, no help available as yet\n");
	exit(0);
}

float gc_content (void) {
	int i;
	int count[5];
	
	for (i = 0; i < 5; i++) count[i] = 0;
	for (i = 0; i < DNA->length; i++) count[(int)DNA->s5[i]]++;
	return (float)(count[1]+count[2]) / (float)(count[0]+count[1]+count[2]+count[3]);
	
}

float repeat_content (void) {
	int i;
	int count = 0;
	for (i = 0; i < DNA->length; i++) {
		switch (DNA->seq[i]) {
			case 'N': case 'n': case 'a': case 'c': case 'g': case 't': count++;
		}
	}
	return (float)count / (float)DNA->length;
}

/*****************************************************************************\
 generic data processing functions
\*****************************************************************************/

void openData (void) {
	DNA_stream = zoeOpenFile(DNA_file);
	ANN_stream = zoeOpenFile(ANN_file);
}

void closeData (void) {
	int    label;
	
	zoeCloseFile(DNA_stream);
	zoeCloseFile(ANN_stream);

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
	for (label = 0; label < zoeLABELS; label++) {
		if (Features[label]) {
			zoeDeleteFeatureVec(Features[label]);
			Features[label] = NULL;
		}
	}
	if (GENE) {
		zoeDeleteCDS(GENE);
		GENE = NULL;
	}
}

int getData (void) {
	zoeFastaFile ff;
	zoeVec genes;
	int i;
	
	/* read new dna and ann */
	if ((ff = zoeReadFastaFile(DNA_stream.stream)) == NULL) return 0;
	
	/* delete old dna and ann objects */
	if (DNA)  zoeDeleteDNA(DNA);
	if (ANTI) zoeDeleteDNA(ANTI);
	if (ANN)  zoeDeleteFeatureTable(ANN);
	if (GENE) zoeDeleteCDS(GENE);
	for (i = 0; i < zoeLABELS; i++) if (Features[i]) zoeDeleteFeatureVec(Features[i]);
	
	/* create new objects */
	DNA = zoeNewDNA(ff->def, ff->seq);
	zoeDeleteFastaFile(ff);
	if (zoeOption("-lcmask")) zoeLCsmooth(DNA, 10, 10, 100);
	ANTI = zoeAntiDNA(DNA->def, DNA);
	ANN  = zoeReadFeatureTable(ANN_stream.stream);

	/* process annotation data */	
	genes = zoeGetGenes(ANN, DNA);
	if (genes->size != 1) {
		zoeExit("must be one GENE per sequence");
	}
	GENE = genes->elem[0];
	zoeDeleteVec(genes);
	
	Features[Acceptor] = extractAcceptor();
	Features[Donor]    = extractDonor();
	Features[Start]    = extractStart();
	Features[Stop]     = extractStop();
	Features[Coding]   = extractCoding();
	Features[Intron]   = extractIntron();
	Features[Inter]    = extractInter();
	Features[UTR5]     = extractUTR5();
	Features[UTR3]     = extractUTR3();
	Features[PolyA]    = extractPolyA();
		
	return 1;
}

zoeVec define_models (void) {
	int i, j;
	struct model * m;
	zoeVec         vec = zoeNewVec();
	float pseudo = 1;
	
	int acc[4]    = {15, 20, 30, 40};
	int don[4]    = { 9, 12, 15, 20};
	int start[4]  = {12, 15, 18, 21};
	int stop[4]   = {9, 12, 15, 18};
	int intron[5] = {1, 2, 3, 4, 5};
	int inter[5]  = {1, 2, 3, 4, 5};
	int coding[4] = {2, 3, 4, 5};
	int utr5[4]   = {1, 2, 3, 4};
	int utr3[4]   = {1, 2, 3, 4};
	int polyA[4]  = {20, 30, 40, 50};
	
	if (VERBOSE) zoeE("defining models:");
	
	/* pseudocounts */
	if (zoeOption("-pseudocount"))  pseudo = atof(zoeOption("-pseudocount"));
	
	/* acceptor */
	if (VERBOSE) zoeE("a");
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 3; j++) {
			m = zoeMalloc(sizeof(struct model));
			m->label  = Acceptor;
			m->length = acc[i];
			m->order  = j;
			m->counts = zoeNewAcceptorModel(j, acc[i], pseudo);
			m->scores = zoeNewAcceptorModel(j, acc[i], 0);
			zoePushVec(vec, m);
		}
	}
	
	/* donor */
	if (VERBOSE) zoeE("d");
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 3; j++) {
			m = zoeMalloc(sizeof(struct model));
			m->label  = Donor;
			m->length = don[i];
			m->order  = j;
			m->counts = zoeNewDonorModel(j, don[i], pseudo);
			m->scores = zoeNewDonorModel(j, don[i], 0);
			zoePushVec(vec, m);
		}
	}
	
	/* start */
	if (VERBOSE) zoeE("m");
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 3; j++) {
			m = zoeMalloc(sizeof(struct model));
			m->label  = Start;
			m->length = start[i];
			m->order  = j;
			m->counts = zoeNewStartModel(j, start[i], pseudo);
			m->scores = zoeNewStartModel(j, start[i], 0);
			zoePushVec(vec, m);
		}
	}
	
	/* stop */
	if (VERBOSE) zoeE("s");
	for (i = 0; i < 4; i++) {
		m = zoeMalloc(sizeof(struct model));
		m->label  = Stop;
		m->length = stop[i];
		m->order  = 0;
		m->counts = zoeNewStopModel(stop[i], pseudo);
		m->scores = zoeNewStopModel(stop[i], 0);
		zoePushVec(vec, m);
	}
	
	/* intron */
	if (VERBOSE) zoeE("i");
	for (i = 0; i < 5; i++) {
		m = zoeMalloc(sizeof(struct model));
		m->label  = Intron;
		m->length = intron[i] +1;
		m->order  = intron[i];
		m->counts = zoeNewIntronModel(intron[i], pseudo);
		m->scores = zoeNewIntronModel(intron[i], 0);
		zoePushVec(vec, m);
	}
	
	/* inter */
	if (VERBOSE) zoeE("n");
	for (i = 0; i < 5; i++) {
		m = zoeMalloc(sizeof(struct model));
		m->label  = Inter;
		m->length = inter[i] +1;
		m->order  = inter[i];
		m->counts = zoeNewInterModel(inter[i], pseudo);
		m->scores = zoeNewInterModel(inter[i], 0);
		zoePushVec(vec, m);
	}
	
	/* coding */
	if (VERBOSE) zoeE("c");
	for (i = 0; i < 4; i++) {
		m = zoeMalloc(sizeof(struct model));
		m->label  = Coding;
		m->length = coding[i] +1;
		m->order  = coding[i];
		m->counts = zoeNewCodingModel(coding[i], pseudo);
		m->scores = zoeNewCodingModel(coding[i], 0);
		
		zoePushVec(vec, m);
	}
	
	/* UTR5 */
	if (VERBOSE) zoeE("5");
	for (i = 0; i < 4; i++) {
		m = zoeMalloc(sizeof(struct model));
		m->label  = UTR5;
		m->length = utr5[i] +1;
		m->order  = utr5[i];
		m->counts = zoeNewUTR5Model(utr5[i], pseudo);
		m->scores = zoeNewUTR5Model(utr5[i], 0);
		zoePushVec(vec, m);
	}
	
	/* UTR3 */
	if (VERBOSE) zoeE("3");
	for (i = 0; i < 4; i++) {
		m = zoeMalloc(sizeof(struct model));
		m->label  = UTR3;
		m->length = utr3[i] +1;
		m->order  = utr3[i];
		m->counts = zoeNewUTR3Model(utr3[i], pseudo);
		m->scores = zoeNewUTR3Model(utr3[i], 0);
		zoePushVec(vec, m);
	}
	
	/* PolyA */
	if (VERBOSE) zoeE("p");
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 3; j++) {
			m = zoeMalloc(sizeof(struct model));
			m->label  = PolyA;
			m->length = polyA[i];
			m->order  = j;
			m->counts = zoeNewPolyAModel(j, polyA[i], pseudo);
			m->scores = zoeNewPolyAModel(j, polyA[i], 0);
			zoePushVec(vec, m);
		}
	}
	
	/* ambiguate counts - deleted later */
	if (VERBOSE) zoeE("...");
	for (i = 0; i < vec->size; i++) {
		m = vec->elem[i];
		zoeAmbiguateModel(m->counts, 0);
	}
	
	if (VERBOSE) zoeE("\n");
	
	return vec;
}



