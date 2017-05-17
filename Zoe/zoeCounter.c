/******************************************************************************\
 zoeCounter - part of the ZOE library for genomic analysis
 
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

\******************************************************************************/

#ifndef ZOE_COUNTER_C
#define ZOE_COUNTER_C

#include "zoeCounter.h"

/******************************************************************************\
 PRIVATE FUNCTIONS
\******************************************************************************/

static char seq2sig (char c) {
	switch (c) {
		case 'A': case 'a':	return  8; /* 1000 */
		case 'C': case 'c':	return  4; /* 0100 */
		case 'G': case 'g':	return  2; /* 0010 */
		case 'T': case 't':	return  1; /* 0001 */
		case 'R': case 'r':	return 10; /* 1010 */
		case 'Y': case 'y':	return  5; /* 0101 */
		case 'M': case 'm':	return 12; /* 1100 */
		case 'K': case 'k':	return  3; /* 0011 */
		case 'W': case 'w':	return  9; /* 1001 */
		case 'S': case 's':	return  6; /* 0110 */
		case 'B': case 'b':	return  7; /* 0111 */
		case 'D': case 'd':	return 11; /* 1011 */
		case 'H': case 'h':	return 13; /* 1101 */
		case 'V': case 'v':	return 14; /* 1110 */
		case 'N': case 'n':	return 15; /* 1111 */
		default:
			zoeWarn("illegal symbols seq2sig");
			return 15;
	}
}

static int zoeSDTlookup (zoeCounter counter, coor_t mfocus) {
	char c;
	int i, j, found, counter_number = -1; /* impossible value */
	
	if (mfocus >= 0) {
		for (i = 0; i < counter->model->submodels; i++) {
			found = 1;
			for (j = 0; j < counter->model->length; j++) {
				if (counter->subcounter[i]->sig[j] == 15) continue;
				c = counter->dna->s16[j + mfocus] | counter->subcounter[i]->sig[j];
				if (c != counter->subcounter[i]->sig[j]) {
					found = 0;
					break;
				}
			}
			if (found) {counter_number = i; break;}
		}
		
	} else {
		/* reposition mfocus */
		mfocus = counter->dna->length -1 + mfocus;
		for (i = 0; i < counter->model->submodels; i++) {
			found = 1;
			for (j = 0; j < counter->model->length; j++) {
				if (counter->subcounter[i]->sig[j] == 15) continue;
				/* use counter->anti rather than counter->dna */
				c = counter->anti->s16[j + mfocus] | counter->subcounter[i]->sig[j];
				if (c != counter->subcounter[i]->sig[j]) {
					found = 0;
					break;
				}
			}
			if (found) {counter_number = i; break;}
		}
	}

	if (counter_number == -1) zoeExit("no counter found? zoeCountSDT");
	return counter_number;
}

static void zoeCountWMM(zoeCounter counter, coor_t pos) {
	coor_t i, mfocus, index;
		
	if (pos >= 0) {
		if ((pos < counter->min_pos) || (pos > counter->max_pos)) return;
		mfocus = pos - counter->model->focus;
		for (i = 0; i < counter->model->length; i++) {
			index = (i * counter->model->symbols) + counter->dna->s5[i + mfocus];
			counter->model->data[index]++;
		}
	} else {
		if ((-pos < counter->min_pos) || (-pos > counter->max_pos)) return;
		/* reposition mfocus */
		mfocus = counter->dna->length -1 + pos - counter->model->focus;
		for (i = 0; i < counter->model->length; i++) {
			/* use anti instead of dna */
			index = (i * counter->model->symbols) + counter->anti->s5[i + mfocus];
			counter->model->data[index]++;
		}
	}
}

static void zoeCountLUT(zoeCounter counter, coor_t pos) {
	coor_t i, p, index, mfocus;
		
	index = 0;
	if (pos >= 0) {
		if ((pos < counter->min_pos) || (pos > counter->max_pos)) return;
		mfocus = pos - counter->model->focus;
		for (i = 0; i < counter->model->length; i++) {
			p = zoePOWER[counter->model->symbols][counter->model->length -i -1];
			index += (p * counter->dna->s5[i + mfocus]);
		}
	} else {
		if ((-pos < counter->min_pos) || (-pos > counter->max_pos)) return;
		/* reposition mfocus */
		mfocus = counter->dna->length -1 + pos - counter->model->focus;
		for (i = 0; i < counter->model->length; i++) {
			p = zoePOWER[counter->model->symbols][counter->model->length -i -1];
			/* use anti instead of dna */
			index += (p * counter->anti->s5[i + mfocus]);
		}
	}
	counter->model->data[index]++;
}

static void zoeCountSAM(zoeCounter counter, coor_t pos) {
	coor_t i, mfocus;
		
	if (pos >= 0) {
		if ((pos < counter->min_pos) || (pos > counter->max_pos)) return;
		mfocus = pos - counter->model->focus;
		
		for (i = 0; i < counter->model->length; i++) {
			counter->subcounter[i]->count(counter->subcounter[i], i + mfocus);
		}
	} else {
		pos = -pos;
		if ((pos < counter->min_pos) || (pos > counter->max_pos)) return;
		
		mfocus = counter->dna->length -1 + pos - counter->model->focus;
		for (i = 0; i < counter->model->length; i++) {
			counter->subcounter[i]->count(counter->subcounter[i], -i - mfocus);
		}
	}
	
}

static void zoeCountSDT(zoeCounter counter, coor_t pos) {
	int counter_number, mfocus;
		
	if (pos >= 0) {if (( pos < counter->min_pos) || ( pos > counter->max_pos)) return;}
	else          {if ((-pos < counter->min_pos) || (-pos > counter->max_pos)) return;}
		
	mfocus = pos - counter->model->focus;
	counter_number = zoeSDTlookup(counter, mfocus);
	counter->subcounter[counter_number]->count(counter->subcounter[counter_number], pos);
}

static void zoeCountMIX(zoeCounter counter, coor_t pos) {
	zoeExit("zoeCountMIX is not yet enabled, need mix ratio");
}

static void zoeCountTRM(zoeCounter counter, coor_t pos) {
	return;
}

static void zoeCountFeature (zoeCounter counter, zoeFeature f) {
	coor_t i;
	
	if (f->strand == '+') {
		for (i = f->start; i <= f->end; i++) {
			counter->count(counter, i);
		}
	} else {
		for (i = f->start; i <= f->end; i++) {
			counter->count(counter, -i);
		}
	}
}

static void zoeCountCDS (zoeCounter counter, zoeFeature f) {
	coor_t  i;
	int     n; /* should be frame_t but screw casting */
	int     start = 0, end = 0;
	
	if (f->strand == '+') {
		switch (f->label) {
			case Einit: case Exon:
				start = f->start +6; /* make room for 5th order Markov Model */
				end   = f->end;      /* this is a hard-coded hack */
				break;
			case Eterm: case Esngl:
				start = f->start +6;
				end   = f->end   -3; /* don't include stop codon */
				break;
			default: zoeExit("attempt to score CDS with non-exon");
		}
		for (i = start; i <= end; i++) {
			n = (i -start +4 - f->inc5) % 3;		
			counter->subcounter[n]->count(counter->subcounter[n], i);
		}
	} else {
		switch (f->label) {
			case Einit: case Exon:
				start = f->start;    /* make room for 5th order Markov Model */
				end   = f->end   -6; /* this is a hard-coded hack */
				break;
			case Eterm: case Esngl:
				start = f->start +3;/* don't include stop codon */
				end   = f->end   -6; 
				break;
			default: zoeExit("zoeCountCDS: attempt to count CDS with non-exon");
		}
		for (i = start; i <= end; i++) {
			n = (i -start +3 - f->inc5) % 3;
			switch (n) {
				case 0: n = 0; break;
				case 1: n = 2; break;
				case 2: n = 1; break;
			}
			counter->subcounter[n]->count(counter->subcounter[n], -i); /* use neg coor */
		}
	}

}

static void zoeIllegalCount (zoeCounter s, coor_t p) {
	zoeExit("illegal count (%s)", s->model->name);
}

/******************************************************************************\
 PUBLIC FUNCTIONS
\******************************************************************************/

void zoeDeleteCounter(zoeCounter counter) {
	int i;
	
	if (counter == NULL) return;
	
	if (counter->model->submodels) {
		for (i = 0; i < counter->model->submodels; i++) {
			zoeDeleteCounter(counter->subcounter[i]);
			counter->subcounter[i] = NULL;
		}
		zoeFree(counter->subcounter);
		counter->subcounter = NULL;
	}
	if (counter->sig) {
		zoeFree(counter->sig);
		counter->sig = NULL;
	}
	counter->model = NULL;
	counter->dna   = NULL;
	counter->anti  = NULL;
	
	zoeFree(counter);
	counter        = NULL;
}

zoeCounter zoeNewCounter(zoeDNA dna, zoeDNA anti, zoeModel model) {
	int        i, j;
	zoeCounter counter = zoeMalloc(sizeof(struct zoeCounter));

	/* determine scoring region */
	counter->min_pos = model->length;
	counter->max_pos = dna->length - model->length -1;
	
	/* set dna and model, clear pointers */
	counter->dna         = dna;
	counter->anti        = anti;
	counter->model       = model;
	counter->subcounter  = NULL;
	counter->sig         = NULL;
	counter->count       = NULL;
	counter->countf      = NULL;
	
	/* bind scoring and counting functions to type of model */
	switch (model->type) {
		case WMM: counter->count = zoeCountWMM; break;
		case LUT: counter->count = zoeCountLUT; break;
		case SAM: counter->count = zoeCountSAM; break;
		case SDT: counter->count = zoeCountSDT; break;
		case MIX: counter->count = zoeCountMIX; break;
		case TRM: counter->count = zoeCountTRM; break;
		default:  counter->count = zoeIllegalCount;
	}
	
	/* bind range counting and scoring functions */
	switch (model->type) {
		case CDS: counter->countf = zoeCountCDS; break;
		default:  counter->countf = zoeCountFeature;		
	}
	
	/* create subcounters for constructed types */
	if (model->type == WMM || model->type == LUT || model->type == TRM) {
		counter->subcounter = NULL;
	} else {
		counter->subcounter = zoeMalloc(model->submodels * sizeof(zoeCounter));
		for (i = 0; i< model->submodels; i++)
			counter->subcounter[i] = zoeNewCounter(dna, anti, model->submodel[i]);
	}
	
	/* create decision tree for SDT */
	if (model->type == SDT) {
		
		/* create sigs for submodels */
		counter->sig = zoeMalloc(model->submodels);
		for (i = 0; i < model->submodels; i++) {
			counter->subcounter[i]->sig = zoeMalloc(model->length);
			for (j = 0; j < model->length; j++) {
				counter->subcounter[i]->sig[j] =
					seq2sig(model->submodel[i]->name[j]);
			}
		}
			
	} else {
		counter->sig = NULL;
	}
	
	return counter;
}


#endif

