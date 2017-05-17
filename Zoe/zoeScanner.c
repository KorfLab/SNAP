/******************************************************************************\
 zoeScanner.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_SCANNER_C
#define ZOE_SCANNER_C

#include "zoeScanner.h"

/******************************************************************************\
 PRIVATE FUNCTIONS
\******************************************************************************/

static char seq2sig (int c) {
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

static int zoeSDTlookup (zoeScanner scanner, coor_t mfocus) {
	char c;
	int i, j, found, scanner_number = -1; /* impossible value */
	
	if (mfocus >= 0) {
		for (i = 0; i < scanner->model->submodels; i++) {
			found = 1;
			for (j = 0; j < scanner->model->length; j++) {
				if (scanner->subscanner[i]->sig[j] == 15) continue;
				c = scanner->dna->s16[j + mfocus] | scanner->subscanner[i]->sig[j];
				if (c != scanner->subscanner[i]->sig[j]) {
					found = 0;
					break;
				}
			}
			if (found) {scanner_number = i; break;}
		}
		
	} else {
		/* reposition mfocus */
		mfocus = scanner->dna->length -1 + mfocus;
		for (i = 0; i < scanner->model->submodels; i++) {
			found = 1;
			for (j = 0; j < scanner->model->length; j++) {
				if (scanner->subscanner[i]->sig[j] == 15) continue;
				/* use scanner->anti rather than scanner->dna */
				c = scanner->anti->s16[j + mfocus] | scanner->subscanner[i]->sig[j];
				if (c != scanner->subscanner[i]->sig[j]) {
					found = 0;
					break;
				}
			}
			if (found) {scanner_number = i; break;}
		}
	}

	if (scanner_number == -1) zoeExit("no scanner found? zoeCountSDT");
	return scanner_number;
}

static score_t zoeScoreWMM (const zoeScanner scanner, coor_t pos) {
	coor_t  i, mfocus, index;
	score_t score, s;
	
	/* user defines and boundaries */
	if (pos >= 0) {
		if (scanner->uscore && scanner->uscore[pos] != MIN_SCORE) return scanner->uscore[pos];
		if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) return MIN_SCORE;
	} else {
		if (scanner->ascore && scanner->ascore[-pos] != MIN_SCORE) return scanner->ascore[-pos];
		if ((-pos < scanner->min_pos) || (-pos > scanner->max_pos)) return MIN_SCORE;
	}
	
	/* scoring */
	score = 0;
	if (pos >= 0) {
		mfocus = pos - scanner->model->focus;
		for (i = 0; i < scanner->model->length; i++) {
			index = (i * scanner->model->symbols) + scanner->dna->s5[i + mfocus];
			s = scanner->model->data[index];
			if (s == MIN_SCORE) return MIN_SCORE;
			score += s;
		}
	} else {
		mfocus = scanner->dna->length -1 + pos - scanner->model->focus;
		for (i = 0; i < scanner->model->length; i++) {
			index = (i * scanner->model->symbols) + scanner->anti->s5[i + mfocus];
			s = scanner->model->data[index];
			if (s == MIN_SCORE) return MIN_SCORE;
			score += s;
		}
	}
	
	/* score adjustment - only used in WMM */
	score += scanner->model->score;
	return score;
}

static score_t zoeScoreLUT (const zoeScanner scanner, coor_t pos) {
	coor_t i, p, index, mfocus;
	
	/* user defines and boundaries */
	if (pos >= 0) {
		if (scanner->uscore && scanner->uscore[pos] != MIN_SCORE) return scanner->uscore[pos];
		if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) return MIN_SCORE;
	} else {
		if (scanner->ascore && scanner->ascore[-pos] != MIN_SCORE) return scanner->ascore[-pos];
		if ((-pos < scanner->min_pos) || (-pos > scanner->max_pos)) return MIN_SCORE;
	}
	
	/* scoring */
	index = 0;
	if (pos >= 0) {
		mfocus = pos - scanner->model->focus;
		for (i = 0; i < scanner->model->length; i++) {
			p = zoePOWER[scanner->model->symbols][scanner->model->length -i -1];
			index += (p * scanner->dna->s5[i + mfocus]);
		}
	} else {
		mfocus = scanner->dna->length -1 + pos - scanner->model->focus;
		for (i = 0; i < scanner->model->length; i++) {
			p = zoePOWER[scanner->model->symbols][scanner->model->length -i -1];
			/* use anti instead of dna */
			index += (p * scanner->anti->s5[i + mfocus]);
		}
	}
	return scanner->model->data[index];	
}

static score_t zoeScoreSAM (const zoeScanner scanner, coor_t pos) {
	coor_t  i;
	score_t score, s;
	coor_t  mfocus;
	
	/* user defines and boundaries */
	if (pos >= 0) {
		if (scanner->uscore && scanner->uscore[pos] != MIN_SCORE) return scanner->uscore[pos];
		if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) return MIN_SCORE;
	} else {
		if (scanner->ascore && scanner->ascore[-pos] != MIN_SCORE) return scanner->ascore[-pos];
		if ((-pos < scanner->min_pos) || (-pos > scanner->max_pos)) return MIN_SCORE;
	}
	
	/* scoring */
	score = 0;
	if (pos >= 0) {
		mfocus = pos - scanner->model->focus;
		for (i = 0; i < scanner->model->length; i++) {
			s = scanner->subscanner[i]->score(scanner->subscanner[i], i + mfocus);
			if (s == MIN_SCORE) return MIN_SCORE;
			score += s;
		}
	} else {
		mfocus = scanner->dna->length -1 + pos - scanner->model->focus;
		for (i = 0; i < scanner->model->length; i++) {
			s = scanner->subscanner[i]->score(scanner->subscanner[i], -i - mfocus);
			if (s == MIN_SCORE) return MIN_SCORE;
			score += s;
		}
	}

	return score;
}

static score_t zoeScoreSDT (const zoeScanner scanner, coor_t pos) {
	int scanner_number, mfocus;
	
	/* user defines and boundaries */
	if (pos >= 0) {
		if (scanner->uscore && scanner->uscore[pos] != MIN_SCORE) return scanner->uscore[pos];
		if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) return MIN_SCORE;
	} else {
		if (scanner->ascore && scanner->ascore[-pos] != MIN_SCORE) return scanner->ascore[-pos];
		if ((-pos < scanner->min_pos) || (-pos > scanner->max_pos)) return MIN_SCORE;
	}
	
	mfocus = pos - scanner->model->focus;
	scanner_number = zoeSDTlookup(scanner, mfocus);
	return scanner->subscanner[scanner_number]->score(scanner->subscanner[scanner_number], pos);
}

static score_t zoeScoreMIX (const zoeScanner scanner, coor_t pos) {
	zoeExit("zoeScoreMix is not yet enabled");
	return 0;
}

static score_t zoeScoreTRM (const zoeScanner scanner, coor_t pos) {
	return MIN_SCORE;
}

static score_t zoeScoreFeature (const zoeScanner scanner, zoeFeature f) {
	coor_t  i;
	score_t s, score = 0;
	
	if (f->strand == '+') {
		for (i = f->start; i <= f->end; i++) {
			s = scanner->score(scanner, i);
			if (s == MIN_SCORE) return MIN_SCORE; /* boundary condition - changed continue */
			score += s;
		}
	} else {
		for (i = f->start; i <= f->end; i++) {
			s = scanner->score(scanner, -i);
			if (s == MIN_SCORE) return MIN_SCORE;
			score += s;
		}
	}
	
	return score;
}

static score_t zoeScoreCDS (const zoeScanner scanner, zoeFeature f) {
	coor_t  i;
	score_t score, s;
	int     n; /* should be frame_t but screw casting */
	int     start = -1, end = -1;
	/*
	zoeDNA     dna = NULL;
	zoeProtein pro;
	*/
	
	if (f->strand == '+') {
		switch (f->label) {
			case Einit: case Exon:
				start = f->start +6; /* make room for 5th order Markov Model */
				end   = f->end;
				break;
			case Eterm: case Esngl:
				start = f->start +6;
				end   = f->end   -3; /* don't include stop codon */
				break;
			case Coding:
				start = f->start +6;
				end   = f->end   -3; /* don't include stop codon (can't be sure) */
				break;
			default:
				zoeExit("zoeScoreCDS: attempt to score CDS with non-coding");
		}
	
		score = 0;
		for (i = start; i <= end; i++) {
			n = (i - start + 4 - f->inc5) % 3;
			s = scanner->subscanner[n]->score(scanner->subscanner[n], i);
			if (s == MIN_SCORE) {
				/*score = MIN_SCORE;*/
				score = 0; /* bizarre, rare error... */
				break;
			}
			
			score += s;
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
			case Coding:
				start = f->start +3; /* don't include stop codon (can't be sure) */
				end   = f->end   -6;
				break;
			default:
				zoeExit("zoeScoreCDS: attempt to score CDS with non-coding");
		}
	
		score = 0;
		for (i = start; i <= end; i++) {
			n = (i -start +3 - f->inc5) % 3;
			switch (n) {
				case 0: n = 0; break;
				case 1: n = 2; break;
				case 2: n = 1; break;
			}
			s = scanner->subscanner[n]->score(scanner->subscanner[n], -i);
			if (s == MIN_SCORE) {
				/*score = MIN_SCORE;*/
				score = 0; /* bizarre, rare error... */
				break;
			}
			score += s;
		}
	}
	
	if (score == MIN_SCORE) {
		zoeWriteFeature(stdout, f);
		
		zoeE("\n\nDarn! Unfortunately there is a bug in the program.\n");
		zoeE("Please send a report to Ian Korf (iankorf@mac.com).\n");
		exit(1);
	}
		
	return score;
}

static score_t zoeIllegalScore (const zoeScanner s, coor_t p) {
	zoeExit("illegal score (%s)", s->model->name);
	return 0;
}

/******************************************************************************\
 PUBLIC FUNCTIONS
\******************************************************************************/

void zoeDeleteScanner(zoeScanner scanner) {
	int i;
	
	if (scanner == NULL) return;
	
	if (scanner->model->submodels) {
		for (i = 0; i < scanner->model->submodels; i++) {
			if (scanner->subscanner[i]) {
				zoeDeleteScanner(scanner->subscanner[i]);
				scanner->subscanner[i] = NULL;
			}
		}
		if (scanner->subscanner) {
			zoeFree(scanner->subscanner);
			scanner->subscanner = NULL;
		}
	}
	if (scanner->sig) {
		zoeFree(scanner->sig);
		scanner->sig = NULL;
	}
	if (scanner->uscore) {
		zoeFree(scanner->uscore);
		scanner->uscore = NULL;
	}
	if (scanner->ascore) {
		zoeFree(scanner->ascore);
		scanner->ascore = NULL;
	}
	
	scanner->model = NULL;
	scanner->dna   = NULL;
	scanner->anti  = NULL;
	
	zoeFree(scanner);
	scanner        = NULL;
}

zoeScanner zoeNewScanner(zoeDNA dna, zoeDNA anti, zoeModel model) {
	int        i, j;
	zoeScanner scanner = zoeMalloc(sizeof(struct zoeScanner));

	/* determine scoring region */
	scanner->min_pos = model->length;
	scanner->max_pos = dna->length - model->length -1;
	
	/* set dna and model, clear pointers */
	scanner->dna         = dna;
	scanner->anti        = anti;
	scanner->model       = model;
	scanner->subscanner  = NULL;
	scanner->sig         = NULL;
	scanner->uscore      = NULL;
	scanner->ascore      = NULL;
	scanner->score       = NULL;
	scanner->scoref      = NULL;
	
	/* bind scoring and counting functions to type of model */
	switch (model->type) {
		case WMM: scanner->score = zoeScoreWMM; break;
		case LUT: scanner->score = zoeScoreLUT; break;
		case SAM: scanner->score = zoeScoreSAM; break;
		case SDT: scanner->score = zoeScoreSDT; break;
		case MIX: scanner->score = zoeScoreMIX; break;
		case TRM: scanner->score = zoeScoreTRM; break;
		default:  scanner->score = zoeIllegalScore;
	}
	
	/* bind range scoring functions */
	switch (model->type) {
		case CDS: scanner->scoref = zoeScoreCDS; break;
		default:  scanner->scoref = zoeScoreFeature;
	}
	
	/* create subscanners for constructed types */
	if (model->type == WMM || model->type == LUT || model->type == TRM) {
		scanner->subscanner = NULL;
	} else {
		scanner->subscanner = zoeMalloc(model->submodels * sizeof(struct zoeScanner));
		for (i = 0; i< model->submodels; i++)
			scanner->subscanner[i] = zoeNewScanner(dna, anti, model->submodel[i]);
	}
	
	/* ensure CDS model is correctly used */
	if (model->type == CDS) {
		if (model->submodels != 3) zoeExit("CDS must have 3 submodels");
		for (i = 0; i < model->submodels; i++) {
			if (model->submodel[i]->type != LUT) zoeExit("CDS submodels not LUT");
		}
	}
	
	/* create decision tree for SDT */
	if (model->type == SDT) {
		
		/* create sigs for submodels */
		scanner->sig = NULL;
		for (i = 0; i < model->submodels; i++) {
			scanner->subscanner[i]->sig = zoeMalloc(model->length);
			for (j = 0; j < model->length; j++) {
				scanner->subscanner[i]->sig[j] =
					seq2sig(model->submodel[i]->name[j]);
			}
		}
			
	} else {
		scanner->sig = NULL;
	}
	
	return scanner;
}

void zoeSetScannerScore(zoeScanner scanner, coor_t pos, score_t score) {
	coor_t i;
		
	if (pos >= 0) {
		if (pos >= scanner->dna->length) {
			zoeWarn("zoeSetScannerScore position (%d) out of range", pos);
			return;
		}
		if (!scanner->uscore) {
			scanner->uscore = zoeMalloc(scanner->dna->length * sizeof(score_t));
			for (i = 0; i < scanner->dna->length; i++) {
				scanner->uscore[i] = MIN_SCORE;
			}
		}
		scanner->uscore[pos] = score;
	} else {
		if (-pos >= scanner->anti->length) {
			zoeWarn("zoeSetScannerScore position (%d) out of range", pos);
			return;
		}
		if (!scanner->ascore) {
			scanner->ascore = zoeMalloc(scanner->anti->length * sizeof(score_t));
			for (i = 0; i < scanner->anti->length; i++) {
				scanner->ascore[i] = MIN_SCORE;
			}
		}
		scanner->ascore[-pos] = score;
	}
}

#endif
