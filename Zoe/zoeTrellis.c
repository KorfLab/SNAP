/******************************************************************************\
 zoeTrellis.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_TRELLIS_C
#define ZOE_TRELLIS_C

#include "zoeTrellis.h"

/****************************************************************************\
 PRIVATE FUNCTIONS
\****************************************************************************/

/*
void zoePrintTerminalBase (int num) {
	
	switch (num) {
		case 0: printf("none"); break;
		case 1: printf("A"); break;
		case 2: printf("C"); break;
		case 3: printf("G"); break;
		case 4: printf("T"); break;
		case 5: printf("N"); break;
		
		case 6: printf("AA"); break;
		case 7: printf("AC"); break;
		case 8: printf("AG"); break;
		case 9: printf("AT"); break;
		case 10: printf("AN"); break;
		
		case 11: printf("CA"); break;
		case 12: printf("CC"); break;
		case 13: printf("CG"); break;
		case 14: printf("CT"); break;
		case 15: printf("CN"); break;
		
		case 16: printf("GA"); break;
		case 17: printf("GC"); break;
		case 18: printf("GG"); break;
		case 19: printf("GT"); break;
		case 20: printf("GN"); break;
		
		case 21: printf("TA"); break;
		case 22: printf("TC"); break;
		case 23: printf("TG"); break;
		case 24: printf("TT"); break;
		case 25: printf("TN"); break;
		
		case 26: printf("NA"); break;
		case 27: printf("NC"); break;
		case 28: printf("NG"); break;
		case 29: printf("NT"); break;
		case 30: printf("NN"); break;
		
		default: zoeExit("impossible error in zoePrintTerminalBase");
	}
}
*/

static int PADDING = 48;

static int PROGRESS_METER = 1;

struct my_max {
	zoeLabel state;
	coor_t   coor;
	score_t  score;
};

static void delete_external_features (zoeTrellis trellis) {
	int state;
	for (state = 0; state < zoeLABELS; state++) {
		if (trellis->features[state] == NULL) continue;
		zoeDeleteFeatureVec(trellis->features[state]);
		trellis->features[state] = NULL;
	}
}

static void transfer_exons (zoeTrellis trellis, zoeFeatureVec sfv) {
	int i;
	coor_t length;
	zoeFeature exon;
	
	for (i = 0; i < sfv->size; i++) {
		exon = sfv->elem[i];
		
		
		
		/* min length filter */
		length = exon->end - exon->start + 1;
		if (length < trellis->min_len[exon->label]) continue;
		
		/* padding filter */
		if (exon->start < PADDING) continue;
		if (exon->end >= trellis->dna->length) continue;
		
		/* create if necessary */
		if (trellis->features[exon->label] == NULL) {
			trellis->features[exon->label] = zoeNewFeatureVec();
		}
		
		zoePushFeatureVec(trellis->features[exon->label], exon);
	}
}

static void compute_external_features (zoeTrellis trellis, coor_t pos) {
	zoeFeatureFactory factory;
	zoeFeatureVec     sfv;
	int               state;
	
	for (state = 0; state < zoeLABELS; state++) {
		if (trellis->factory[state] == NULL) continue;
		factory = trellis->factory[state];
		if (factory == NULL) continue;

		sfv = factory->create(factory, pos);
		if (sfv == NULL) continue;
		
		switch (state) {
			case Exon:
				transfer_exons(trellis, sfv);
				zoeDeleteFeatureVec(sfv);
				break;
			default:
				if (sfv->last->start < PADDING) sfv->last->start = PADDING;
				trellis->features[state] = sfv;			
		}
	}
}

static int legal_first_jump (zoeLabel int_state, zoeFeature exon) {

	switch (int_state) {
		case Int0:
			if (exon->inc5 == 0) return 1;
			else                  return 0;
		case Int1: case Int1T:
			if (exon->inc5 == 2) return 1;
			else                  return 0;
		case Int2: case Int2TA: case Int2TG:
			if (exon->inc5 == 1) return 1;
			else                  return 0;
		default: return 1;
	}
}

static int legal_second_jump (zoeDNA dna, zoeFeature exon, zoeLabel int_state) {
	int terminal_base;
		
	switch (exon->inc3) {
		case 1:  terminal_base =  dna->s5[exon->end] +1; break;
		case 2:  terminal_base = (dna->s5[exon->end -1] +1) * 5
			                   +  dna->s5[exon->end] +1; break;
		default: terminal_base = 0;
	}

	switch (int_state) {
		case Int0: /* must have no terminal bases */
			switch (terminal_base) {
				case 0:  return 1;
				default: return 0;
			}
		case Int1: /* must have one terminal, but not T */
			switch (terminal_base) {
				case 1: case 2: case 3: case 5: return 1;
				default:                        return 0;
			}
		case Int1T: /* must have terminal T */
			switch (terminal_base) {
				case 4:  return 1;
				default: return 0;
			}
		case Int2: /* must have two terminals, but not TA or TG  */
			if (terminal_base > 5
				&& terminal_base != 21
				&& terminal_base != 23) return 1;
			else                        return 0;
		case Int2TA: /* must have terminal TA */
			switch (terminal_base) {
				case 21: return 1;
				default: return 0;
			}
		case Int2TG: /* must have terminal TG */
			switch (terminal_base) {
				case 23: return 1;
				default: return 0;
			}
		default: return 1;
	}
}

static int creates_stop_codon (zoeDNA dna, zoeLabel pre_state, zoeFeature exon) {
	int initial_base;
	
	switch (exon->inc5) {
		case 1:  initial_base =  dna->s5[exon->start]    + 1; break;
		case 2:  initial_base = (dna->s5[exon->start]    + 1) * 5
		                      +  dna->s5[exon->start +1] + 1; break;
		default: initial_base = 0;
	}
	
	switch (pre_state) {
		case Int1T:
			switch (initial_base) {
				case 6: case 8: case 16: return 1; /* T + {AA, AG, GA} */
				default:                 return 0;
			}
		case Int2TA:
			switch (initial_base) {
				case 1: case 3: return 1; /* TA + {A, G} */
				default:        return 0;
			}
		case Int2TG:
			switch (initial_base) {
				case 1:  return 1; /* TG + A */
				default: return 0;
			}
		default: return 0;
	}
}

static score_t pre_state_duration_score (zoeTrellis trellis, zoeLabel state, int int_end) {
	int i, stop, length;
	
	/* geometric */
	if (trellis->hmm->smap[state]->geometric) /* length 1 geometric cost */
		return zoeScoreDuration(trellis->hmm->dmap[state], 1);
	
	
	/* explicit */
	stop = PADDING;
	for (i = int_end; i > stop; i--) if (trellis->trace[state][i] != -1) break;
	length = int_end -i + 1 + trellis->hmm->cmap[state];
		
	if (length < trellis->min_len[state]) return MIN_SCORE;
	else return zoeScoreDuration(trellis->hmm->dmap[state], length);
	
	/* error to reach */
	zoeExit("should not get here in pre-state duration score");
	return MIN_SCORE;
}

static score_t internal_score (zoeTrellis trellis, zoeLabel state, coor_t i) {
	zoeScanner scanner;
	score_t    content_score, extend_score, prev_score;
	
	if (trellis->score[state][i-1] == MIN_SCORE) return MIN_SCORE; /* keep MIN_SCORE */
	
	scanner       = trellis->scanner[state];
	content_score = scanner->score(scanner, i);
	extend_score  = trellis->hmm->xmap[state];
	prev_score    = trellis->score[state][i-1];

	return content_score + extend_score + prev_score - trellis->exp_score;
}

struct maxExt {
	score_t    score;
	zoeFeature feature;
	zoeLabel   pre_state;
};

static struct maxExt external_score (
	zoeTrellis trellis,
	coor_t     pos,
	zoeLabel   int_state)
{
	int               i, j, length;
	score_t           cscore, dscore, t1score, t2score, phscore1, phscore2,
	                  xscore, total_score, pre_score, exp_score, pro_score = 0;
	zoeLabel          pre_state, ext_state;
	zoeIVec           ivec;
	zoeFeature        f;
	zoeHMM            hmm = trellis->hmm;
	zoeDNA            dna = trellis->dna;
	zoeFeatureVec     sfv;
	struct maxExt     max;
	
	max.score     = MIN_SCORE;
	max.feature   = NULL;
	max.pre_state = -1;
		
/*	    pre_state          ext_state        int_state
	-----------------|||||||||||||||||||||------------
	     intron               exon           intron
	                ->                   ->
	            transition1          transition2
	            
	<---------------|
         xscore
	              ...[......content......]...
	                 <------duration----->
*/

	
	/* external states */	
	for (ext_state = 0; ext_state < zoeLABELS; ext_state++) {
		if (hmm->jmap[int_state][ext_state] == NULL) continue;
						
		ivec = hmm->jmap[int_state][ext_state];
		
		for (i = 0; i < ivec->size; i++) {
			pre_state = ivec->elem[i];
			
			sfv = trellis->features[ext_state];
			if (sfv == NULL) continue;
			
			for (j = 0; j < sfv->size; j++) {
				f = sfv->elem[j];
				
				length    = f->end - f->start +1;
				pre_score = trellis->score[pre_state][pos -length];
				
				/* filters */
				if (pre_score == MIN_SCORE) continue;
				
				if (f->label == Repeat || f->label == ORF || f->label == CNS) {
					if (pre_state != int_state) continue; /* shuttle */
				} else {
					if (!legal_first_jump(pre_state, f))        continue;
					if (!legal_second_jump(dna, f, int_state))  continue;
					if (creates_stop_codon(dna, pre_state, f))  continue;
				}
				
				/* pre-state duration score */
				xscore = pre_state_duration_score(trellis, pre_state, f->start -1);

				if (xscore == MIN_SCORE) continue;
				
				/* expected score */
				exp_score = trellis->exp_score * (f->end - f->start +1);
				f->score -= exp_score;
	
				/* 'exon' scores */
				cscore = f->score;
				dscore = zoeScoreDuration(hmm->dmap[f->label], length);
								
				/* transition scores */
				t1score = hmm->tmap[f->label][int_state];
				t2score = hmm->tmap[pre_state][f->label];
				
				/* phase scores & profile score */
				switch (ext_state) {
					case Einit: case Eterm: case Exon: case Esngl:
						phscore1  = zoeScorePhase(hmm->phasepref, pre_state, ext_state, 0);
						phscore2  = zoeScorePhase(hmm->phasepref, ext_state, int_state, f->inc5);
						if (trellis->ext) {
							pro_score = trellis->ext(trellis, pos, pre_state, f);
							f->score += pro_score;
						}
						break;
					default:
						phscore1  = 0;
						phscore2  = 0;
						pro_score = 0;
				}
				
				/* total score */
				total_score = cscore + dscore + t1score + t2score
					+ xscore + pre_score + phscore1 + phscore2 + pro_score;
							
				if (total_score > max.score) {
					zoeDeleteFeature(max.feature);
					max.score     = total_score;
					max.pre_state = pre_state;
					max.feature   = zoeCopyFeature(f);
				}
			}
		}
	}
	
	return max;
}

static zoeFeatureVec trace_trellis (zoeTrellis trellis, zoeLabel state) {
	zoeFeatureVec sfv;
	int           i, int_end;
	zoeFeature    estate = NULL, istate = NULL;
	int           trace, length;
	score_t       expected;
		
	sfv = zoeNewFeatureVec();
		
	i = trellis->dna->length -1 -PADDING;
	int_end = i;
	while (i > PADDING) {
		i--;
		trace = trellis->trace[state][i];
		if (trace < 0) continue;
		
		/* internal state */
		istate = zoeNewFeature(
			state,
			i,
			int_end,
			'+', 0, 0, 0, 0, NULL/*, NULL*/);
		istate->score = trellis->scanner[state]->scoref(trellis->scanner[state], istate);
		zoePushFeatureVec(sfv, istate);
		zoeDeleteFeature(istate);
		
		/* external state */
		estate = trellis->keep->elem[trace];
		if (estate->label == Repeat && estate->start < PADDING)
			estate->start = PADDING; /* Repeats and PADDING conspire to madness */
		zoePushFeatureVec(sfv, estate);
		
		/* update */
		state   = trellis->jump->elem[trace];
		i       = estate->start -1;
		int_end = i;
	}
	if (sfv->size == 0) return sfv;
	
	istate = zoeNewFeature(
		state, i, int_end, '+',
		trellis->scanner[state]->scoref(trellis->scanner[state], istate),
		0, 0, 0, NULL/*, NULL*/);
	
	zoePushFeatureVec(sfv, istate);
	zoeDeleteFeature(istate);
	
	/* rescore with expected */
	for (i = 0; i < sfv->size; i++) {
		length = sfv->elem[i]->end - sfv->elem[i]->start +1;
		expected = length * trellis->exp_score;
		sfv->elem[i]->score -= expected;
	}
		
	return sfv;
}

static zoeFeatureTable label_genes (const char *def, zoeFeatureVec features) {
	int             i, group;
	char            id[64], name[256];
	zoeFeature      f;
	zoeFeatureTable table = zoeNewFeatureTable(def, NULL);
	
	if (features->size == 0) return table;
	
	strncpy(id, def, 63);
	for (i = 0; i < strlen(id); i++) {
		if (isspace((int)id[i]) || i == 63) {
			id[i] = '\0';
			break;
		}
	}
	
	group = features->last->label == Inter ? 0 : 1;
	sprintf(name, "%s-snap.%d", id, group);
	
	for (i = features->size -1; i >= 0; i--) {
		f = features->elem[i];
						
		switch (f->label) {
			case Inter: case None: /* gene delimiters */
				group++;
				sprintf(name, "%s-snap.%d", id, group);
				zoeAddFeature(table, f);
				break;
			case Esngl: case Eterm:
				f->end += 3;
				/* intentional fall through! */
			case Exon: case Einit:
				f->group = name;
				zoeAddFeature(table, f);
				f->group = NULL;
				break;
			case Int0: case Int1: case Int1T: case Int2: case Int2TA: case Int2TG:
				/* rename all specific introns to general */
				f->label = Intron;
				f->group = name;
				zoeAddFeature(table, f);
				f->group = NULL;
				break;
			default:
				f->group = name;
				zoeAddFeature(table, f);
				f->group = NULL;
				break;
		}
				
	}
		
	return table;
}

static void xdefine_intron (zoeTrellis trellis, const zoeFeature intron) {
	int i;
	score_t s;
	
	if (strcmp(intron->group, "SET") == 0) {
		for (i = intron->start +3; i <= intron->end -3; i++) {
			zoeSetScannerScore(trellis->scanner[Int0],   i, intron->score);
			zoeSetScannerScore(trellis->scanner[Int1],   i, intron->score);
			zoeSetScannerScore(trellis->scanner[Int1T],  i, intron->score);
			zoeSetScannerScore(trellis->scanner[Int2],   i, intron->score);
			zoeSetScannerScore(trellis->scanner[Int2TA], i, intron->score);
			zoeSetScannerScore(trellis->scanner[Int2TG], i, intron->score);
		}
	} else if (strcmp(intron->group, "ADJ") == 0) {
		for (i = intron->start +3; i <= intron->end -3; i++) {
			s = trellis->scanner[Int0]->score(trellis->scanner[Int0], i);
			s += intron->score;
			zoeSetScannerScore(trellis->scanner[Int0],   i, s);
			zoeSetScannerScore(trellis->scanner[Int1],   i, s);
			zoeSetScannerScore(trellis->scanner[Int1T],  i, s);
			zoeSetScannerScore(trellis->scanner[Int2],   i, s);
			zoeSetScannerScore(trellis->scanner[Int2TA], i, s);
			zoeSetScannerScore(trellis->scanner[Int2TG], i, s);
		}
	} else {
		zoeExit("unrecognized command (%s)", intron->group);
	}
}


static void xdefine_coding (zoeTrellis trellis, const zoeFeature coding) {
	int i, j;
	score_t s;
	zoeScanner scan = trellis->scanner[Coding];
	zoeScanner sub;
	
	/*
	
		Note: this function does not discriminate among the various frames. It
		SETs or ADJusts all coding scores in the region regardless of frame.
		Why? Because I'm currently lazy and I haven't the time to debug the
		frame just now. It does not change stop codons from having MIN_SCORE
		however.
	
	*/
	
	
	if (strcmp(coding->group, "SET") == 0) {
		for (i = coding->start +3; i <= coding->end -3; i++) {
			for (j = 0; j < 3; j++) {
				sub = scan->subscanner[j];
				s = sub->score(sub, i);
				if (s != MIN_SCORE) {
					zoeSetScannerScore(sub, i, coding->score);
				}
			}
		}
	} else if (strcmp(coding->group, "ADJ") == 0) {
		for (i = coding->start +3; i <= coding->end -3; i++) {
			for (j = 0; j < 3; j++) {
				sub = scan->subscanner[j];
				s = sub->score(sub, i);
				if (s != MIN_SCORE) {
					zoeSetScannerScore(sub, i, s + coding->score);
				}
			}
		}
	} else {
		zoeExit("unrecognized command (%s)", coding->group);
	}
}


static void adefine_trellis (zoeTrellis trellis, zoeLabel label) {
	int        i, j;
	score_t    score, ascore;
	zoeScanner scanner;
	
	ascore = zoeGetAscore(label);
	
	switch (label) {
		case Acceptor: case Donor: case Start: case Stop: case Inter:
			scanner = trellis->scanner[label];
			for (i = 0; i < trellis->dna->length; i++) {
				score = scanner->score(scanner, i);
				if (score == MIN_SCORE) continue;
				zoeSetScannerScore(scanner, i, score + ascore);
			}
			break;
		case Intron:
			for (j = Int0; j <= Int2TA; j++) {
				scanner = trellis->scanner[j];
				for (i = 0; i < trellis->dna->length; i++) {
					score = scanner->score(scanner, i);
					if (score == MIN_SCORE) continue;
					zoeSetScannerScore(scanner, i, score + ascore);
				}
			}
			break;
		case Coding:
			for (j = 0; j < 3; j++) {
				scanner = trellis->scanner[label]->subscanner[j];
				for (i = 0; i < trellis->dna->length; i++) {
					score = scanner->score(scanner, i);
					if (score == MIN_SCORE) continue;
					zoeSetScannerScore(scanner, i, score + ascore);
				}
			}
			break;
		default: zoeExit("Ascore not yet implemented for %d", label);
	}
	
}

static void xdefine_trellis (zoeTrellis trellis, const zoeFeatureVec vec) {
	int           i, j;
	zoeScanner    scanner;
	zoeFeature    f;
	score_t       score;
	char          name[32];
	
	/* copy all features and pad */
	trellis->xdef = zoeNewFeatureVec();
	for (i = 0; i < vec->size; i++) {
		f = vec->elem[i];
		zoePushFeatureVec(trellis->xdef, f);
		trellis->xdef->last->start += PADDING;
		trellis->xdef->last->end   += PADDING;
	}
	
	/* xdef loop */
	for (i = 0; i < trellis->xdef->size; i++) {
		f = trellis->xdef->elem[i];
		
		/* coordinate checks */
		if (f->start < PADDING || f->end > trellis->dna->length) {
			zoeExit("xdef out of bounds (%d)", f->end -PADDING +1);
		}
		
		/* scanner check */
		scanner = trellis->scanner[f->label];
		if (f->label != Intron && scanner == NULL) {
			zoeLabel2Text(f->label, name);
			zoeWriteFeature(stderr, f);
			zoeExit("attempt to xdef trellis with unknown scanner (%s)", name);
		}
		
		/* command check */
		if (f->group == NULL) {
			zoeWriteFeature(stderr, f);
			zoeExit("xdef has no command");
		}
		
		/* xdef */
		if (f->label == Intron) {
			xdefine_intron(trellis, f);
		} else if (f->label == Coding) {
			xdefine_coding(trellis, f);
		} else if (strcmp(f->group, "SET") == 0) {
			for (j = f->start; j <= f->end; j++) {
				zoeSetScannerScore(scanner, j, f->score);
			}
		} else if (strcmp(f->group, "ADJ") == 0) {
			for (j = f->start; j <= f->end; j++) {
				score = scanner->score(scanner, j);
				score += f->score;
				zoeSetScannerScore(scanner, j, score);
			}
		} else if (strcmp(f->group, "OK") == 0) {
			zoeExit("OK not advised");
			for (j = f->start; j <= f->end; j++) {
				score = scanner->score(scanner, j);
				if (score == MIN_SCORE) {
					zoeSetScannerScore(scanner, j, f->score);
				}
			}
		} else {
			zoeWriteFeature(stderr, f);
			zoeExit("unrecognized xdef command (%s)", f->group);
		}
	}
}

static score_t expected_score (const zoeDNA dna) {
	int     i;
	int     count[5];
	int     total = 0;
	float   freq[4];
	score_t score[4];
	score_t exp_score = 0;
	
	/* init */
	for (i = 0; i < 5; i++) count[i] = 0;
	
	/* nucleotide composition of dna (ignoring Ns) */
	for (i = 0; i < dna->length; i++) count[(int)dna->s5[i]]++;
	total = count[0] + count[1] + count[2] + count[3];
	for (i = 0; i < 4; i++) freq[i] = (float)count[i] / (float)total;
	
	/* expected score */
	for (i = 0; i < 4; i++) score[i] = zoeFloat2Score(freq[i] / 0.25); /* titus */
	for (i = 0; i < 4; i++) exp_score += score[i] * freq[i];
	
	return exp_score;
}

/****************************************************************************\
 PUBLIC FUNCTIONS
\****************************************************************************/


void zoeDeleteTrellis (zoeTrellis trellis) {
	int i;
	
	for (i = 0; i < zoeLABELS; i++) {
		if (trellis->scanner[i]  != NULL) zoeDeleteScanner(trellis->scanner[i]);
		if (trellis->trace[i]    != NULL) zoeFree(trellis->trace[i]);
		if (trellis->score[i]    != NULL) zoeFree(trellis->score[i]);
		if (trellis->factory[i]  != NULL) zoeDeleteFeatureFactory(trellis->factory[i]);
	}
	
	zoeDeleteDNA(trellis->dna);
	zoeDeleteDNA(trellis->anti);
	zoeDeleteFeatureVec(trellis->keep);
	zoeDeleteIVec(trellis->jump);
	zoeFree(trellis);
	
}


zoeTrellis zoeNewTrellis (
	const zoeDNA real_dna,
	const zoeHMM hmm,
	const zoeFeatureVec xdef)
{
	int          i, label;
	zoeState     state;
	zoeDNA       dna, anti;
	zoeTrellis   trellis;
	int          MinimumRepeatLength = 10;
	int          MinimumORFScore = 0;
	
	trellis = zoeMalloc(sizeof(struct zoeTrellis));
		
	/* clear out pointers */
	trellis->dna   = NULL;
	trellis->hmm   = NULL;
	trellis->ext   = NULL;
	for (label = 0; label < zoeLABELS; label++) {
		trellis->scanner[label]  = NULL;
		trellis->trace[label]    = NULL;
		trellis->score[label]    = NULL;
		trellis->factory[label]  = NULL;
		trellis->internal[label] = 0;
		trellis->features[label] = NULL;
	}
	
	/* initial setup */
	dna = zoeMakePaddedDNA(real_dna, PADDING);
	anti = zoeAntiDNA(dna->def, dna);
	trellis->dna       = dna;
	trellis->anti      = anti;
	trellis->hmm       = hmm;
	trellis->xdef      = xdef;
	trellis->max_score = MIN_SCORE;
	trellis->exp_score = expected_score(real_dna);
	trellis->keep      = zoeNewFeatureVec();
	trellis->jump      = zoeNewIVec();

	/* create scanners */
	for (label = 0; label < zoeLABELS; label++) {
		if (hmm->mmap[label] == NULL) continue;
		trellis->scanner[label] = zoeNewScanner(dna, anti, hmm->mmap[label]);
	}
	
	/* modify scanners with arbitrary scoring filters */
	for (label = 0; label < zoeLABELS; label++) {
		if (zoeGetAscore(label) != 0) adefine_trellis(trellis, label);
	}
	
	/* modify scanners with external information */
	if (xdef) xdefine_trellis(trellis, xdef);

	/* create factories for external & shuttle states */
	if (PROGRESS_METER) zoeE("scoring");
	for (i = 0; i <  hmm->states; i++) {
		state = hmm->state[i];
		if (state->type == INTERNAL) continue;
		if (PROGRESS_METER) zoeE(".");
		switch (state->label) {
			case Exon: case Einit: case Eterm: case Esngl:
				if (trellis->factory[Exon]) break; /* set only once */
				trellis->factory[Exon] = zoeNewEFactory(
					trellis->scanner[Coding],
					trellis->scanner[Acceptor],
					trellis->scanner[Donor],
					trellis->scanner[Start],
					trellis->scanner[Stop]);
				break;
			case Repeat:
				trellis->factory[Repeat] = zoeNewRFactory(
					trellis->scanner[Repeat], MinimumRepeatLength);
				break;
			case PolyA:
				trellis->factory[PolyA] = zoeNewSFactory(trellis->scanner[PolyA], PolyA);
				break;
			case Prom:
				trellis->factory[Prom] = zoeNewSFactory(trellis->scanner[Prom], Prom);
				break;
			case TSS:
				trellis->factory[TSS] = zoeNewSFactory(trellis->scanner[TSS], TSS);
				break;
			case ORF:
				trellis->factory[ORF] = zoeNewXFactory(
					trellis->scanner[Coding],
					trellis->scanner[Acceptor],
					trellis->scanner[Donor],
					trellis->scanner[Start],
					trellis->scanner[Stop],
					state->min,
					MinimumORFScore);
				break;
			default:
				zoeExit("zoeNewTrellis: can't make such a factory\n");
				break;
		}
	}
	
	/* trace & scores: 2D matrices */
	for (i = 0; i < hmm->states; i++) {
		state = hmm->state[i];
		if (state->type != INTERNAL) continue;
		
		if (state->label == Intron) {
		
			trellis->internal[Int0]   = 1;
			trellis->internal[Int1]   = 1;
			trellis->internal[Int1T]  = 1;
			trellis->internal[Int2]   = 1;
			trellis->internal[Int2TA] = 1;
			trellis->internal[Int2TG] = 1;
		
			trellis->trace[Int0]   = zoeCalloc(dna->length, sizeof(int));
			trellis->trace[Int1]   = zoeCalloc(dna->length, sizeof(int));
			trellis->trace[Int1T]  = zoeCalloc(dna->length, sizeof(int));
			trellis->trace[Int2]   = zoeCalloc(dna->length, sizeof(int));
			trellis->trace[Int2TA] = zoeCalloc(dna->length, sizeof(int));
			trellis->trace[Int2TG] = zoeCalloc(dna->length, sizeof(int));
			
			trellis->score[Int0]   = zoeCalloc(dna->length, sizeof(score_t));
			trellis->score[Int1]   = zoeCalloc(dna->length, sizeof(score_t));
			trellis->score[Int1T]  = zoeCalloc(dna->length, sizeof(score_t));
			trellis->score[Int2]   = zoeCalloc(dna->length, sizeof(score_t));
			trellis->score[Int2TA] = zoeCalloc(dna->length, sizeof(score_t));
			trellis->score[Int2TG] = zoeCalloc(dna->length, sizeof(score_t));
			
		} else {
		
			trellis->internal[state->label] = 1;
			trellis->trace[state->label]    = zoeCalloc(dna->length, sizeof(int));
			trellis->score[state->label]    = zoeCalloc(dna->length, sizeof(score_t));

		}
	}
	
	/* minimum and maximum lengths */
	for (label = 0; label < zoeLABELS; label++) {
		if (hmm->smap[label] == NULL) continue;
		trellis->min_len[label] = hmm->smap[label]->min;
		trellis->max_len[label] = hmm->smap[label]->max;
	}
	
	return trellis;
}


zoeVec zoePredictGenes (zoeTrellis trellis) {
	coor_t          i;           /* iterator for sequence */
	int             j;           /* iterator for internal states */
	score_t         iscore;      /* score of internal path */
	zoeHMM          hmm = trellis->hmm;
	zoeDNA          dna = trellis->dna;
	zoeFeatureTable table;
	zoeFeatureVec   features;
	zoeVec          genes;
	zoeCDS          gene;
	zoeLabel        max_state;
	score_t         max_score;
	score_t         terminal_score;
	struct maxExt   emax;
	int             progress;
	int             percent;
			
	/*-------------------------------------------------*
	 |             [0] [1] [2] [3] [4] [5] [6] [7] ... |
	 | (None)  [0]                                     |
	 | (Inter) [1]              X                      |
	 | (...)   [2]                                     |
	 | :                  X = trellis->cell[1][3]      |
	 | :                               cell[j][i]      |
	 *-------------------------------------------------*/
	 	 	
	/* initialization */
	for (j = 0; j < zoeLABELS; j++) {
		if (trellis->internal[j] == 0) continue;
		for (i = 0; i <= PADDING; i++) {
			trellis->score[j][i] = hmm->imap[j];
			trellis->trace[j][i]  = -1;
		}
	}
	
	/* induction */
	if (PROGRESS_METER) zoeE("decoding");
	progress = dna->length / 20;
	percent = 0;
	
	for (i = PADDING; i < dna->length - PADDING; i++) {
	
		if (PROGRESS_METER && i % progress == 0) {
			if (i % (progress * 2) == 0) {
				percent += 10;
				zoeE("%d", percent);
			} else {
				zoeE(".");
			}
		}
	
		compute_external_features(trellis, i);
		
		for (j = 0; j < zoeLABELS; j++) {
			if (trellis->internal[j] == 0) continue;
						
			iscore = internal_score(trellis, j, i);
			emax   = external_score(trellis, i, j);
			
			if (iscore == MIN_SCORE && emax.score == MIN_SCORE) {
				trellis->score[j][i] = MIN_SCORE;
				trellis->trace[j][i]  = -1;
			} else if (iscore == MIN_SCORE || emax.score > iscore) {
				trellis->score[j][i] = emax.score;
				trellis->trace[j][i]  = trellis->keep->size;
				zoePushFeatureVec(trellis->keep, emax.feature);
				zoePushIVec(trellis->jump, emax.pre_state);
			} else if (emax.score == MIN_SCORE || emax.score <= iscore) {
				trellis->score[j][i] = iscore;
				trellis->trace[j][i]  = -1;
			} else {
				zoeExit("um, didn't think that was possible");
			}
			
			if (emax.feature) zoeDeleteFeature(emax.feature);
		}
		
		delete_external_features(trellis);
	}
	if (PROGRESS_METER) zoeE("100");
	
	/* find maximum ending state */
	max_state = None;
	max_score = MIN_SCORE;
	for (j = 0; j < zoeLABELS; j++) {
		if (trellis->internal[j] == 0) continue;
		if (hmm->kmap[j] == MIN_SCORE) continue; /* no sense computing */
		terminal_score = hmm->kmap[j] + trellis->score[j][dna->length -1 -PADDING];
		if (terminal_score > max_score) {
			max_state = j;
			max_score = terminal_score;
		}
	}
	if (max_score == MIN_SCORE) zoeExit("traceback from MIN_SCORE");
	trellis->max_score = max_score;
	
	/* get genes from trace-back */
	features = trace_trellis(trellis, max_state);
	table    = label_genes(dna->def, features);
	genes    = zoeGetGenes(table, dna);

			
	/* score all genes - except when using external scoring function */
	if (!trellis->ext) {
		for (i = 0; i < genes->size; i++) {
			gene = genes->elem[i];
			zoeScoreCDS(trellis, gene, 1, 0); /* erasing the external scores! */
		}
	}
	
	/* remove padding */
	for (i = 0; i < genes->size; i++) {
		gene = genes->elem[i];
		gene->start -= PADDING;
		gene->end   -= PADDING;
		for (j = 0; j < gene->exons->size; j++) {
			gene->exons->elem[j]->start -= PADDING;
			gene->exons->elem[j]->end   -= PADDING;
		}
		for (j = 0; j < gene->introns->size; j++) {
			gene->introns->elem[j]->start -= PADDING;
			gene->introns->elem[j]->end   -= PADDING;
		}
		for (j = 0; j < gene->source->size; j++) {
			gene->source->elem[j]->start -= PADDING;
			gene->source->elem[j]->end   -= PADDING;
		}
	}
	
	/* clean up */
	zoeDeleteFeatureVec(features);
	zoeDeleteFeatureTable(table);
	
	if (PROGRESS_METER) zoeE(" done\n");
	return genes;
}

void zoeScoreCDS (zoeTrellis t, zoeCDS cds, int padded, int error_ok) {
	int i;
	
	cds->score = 0;
		
	for (i = 0; i < cds->exons->size; i++) {
		cds->exons->elem[i]->score = zoeScoreExon(t, cds->exons->elem[i], padded, error_ok);
		cds->score += cds->exons->elem[i]->score;
	}
	
	for (i = 0; i < cds->introns->size; i++) {
		cds->introns->elem[i]->score = zoeScoreIntron(t, cds->introns->elem[i], padded);
		cds->score += cds->introns->elem[i]->score;
	}
}

void zoeSetTrellisMeter (int val) {
	PROGRESS_METER = val;
}

void zoeSetTrellisPadding (int val) {
	PADDING = val;
}

char* zoeGetPartialProtein (zoeTrellis trellis, zoeLabel pre_state, zoeFeature last_exon) {
	zoeFeatureVec sfv;
	zoeFeature    feature, exon;
	int           i, j, idx, trace, tx_length;
	char          *tx, *aa;	
	
	/* get exons */
	sfv = zoeNewFeatureVec();
	zoePushFeatureVec(sfv, last_exon);
	i = last_exon->start;
	while (i > PADDING) {
		i--;
		trace = trellis->trace[pre_state][i];
		if (trace < 0) continue;
		feature = trellis->keep->elem[trace];
		switch (feature->label) {
			case Einit:
			case Eterm:
			case Esngl:
			case Exon: zoePushFeatureVec(sfv, feature); break;
			default: break;
		}
		
		pre_state  = trellis->jump->elem[trace];
		i          = feature->start -1;
	}
	if (sfv->size == 0) return NULL;
		
	/* transcribe */
	tx_length = 0;
	for (i = sfv->size -1; i >= 0; i--) {
		exon = sfv->elem[i];
		tx_length += exon->end - exon->start +1;
	}
	tx = zoeMalloc(tx_length+1);
	idx = 0;
	for (i = sfv->size -1; i >= 0; i--) {
		exon = sfv->elem[i];
		for (j = exon->start; j <= exon->end; j++) {
			tx[idx] = trellis->dna->s5[j];
			idx++;
		}
	}
	tx[idx] = '\0';

	/* translate */
	aa = zoeTranslateS5(tx, tx_length, sfv->elem[sfv->size-1]->inc5);
	
	/* clean up */
	zoeDeleteFeatureVec(sfv);
	zoeFree(tx);
	
	return aa;
}

score_t zoeScoreExon (zoeTrellis t, zoeFeature exon, int padded, int error_ok) {
	zoeScanner scan5, scan3;
	score_t    score, cscore, dscore, xscore, score5, score3;

	if (!padded) {
		exon->start += PADDING;
		exon->end   += PADDING;
	}
	
	/* 5'and 3' scores */
	if (exon->strand == '+') {
		switch (exon->label) {
			case Einit:
				scan5 = t->scanner[Start];
				scan3 = t->scanner[Donor];
				score5 = scan5->score(scan5, exon->start);
				score3 = scan3->score(scan3, exon->end +1);
				break;
			case Esngl:
				scan5 = t->scanner[Start];
				scan3 = t->scanner[Stop];
				score5 = scan5->score(scan5, exon->start);
				score3 = scan3->score(scan3, exon->end -2);
				break;
			case Eterm:
				scan5 = t->scanner[Acceptor];
				scan3 = t->scanner[Stop];
				score5 = scan5->score(scan5, exon->start -1);
				score3 = scan3->score(scan3, exon->end -2);
				break;
			case Exon:
				scan5 = t->scanner[Acceptor];
				scan3 = t->scanner[Donor];
				score5 = scan5->score(scan5, exon->start -1);
				score3 = scan3->score(scan3, exon->end +1);
				break;
			default:
				score5 = 0;
				score3 = 0; /* optimizer shush */
				zoeExit("not possible");
		}
	} else {
		switch (exon->label) {
			case Einit:
				scan5 = t->scanner[Start];
				scan3 = t->scanner[Donor];
				score5 = scan5->score(scan5, -exon->end);
				score3 = scan3->score(scan3, -exon->start +1);
				break;
			case Esngl:
				scan5 = t->scanner[Start];
				scan3 = t->scanner[Stop];
				score5 = scan5->score(scan5, -exon->end);
				score3 = scan3->score(scan3, -exon->start -2);
				break;
			case Eterm:
				scan5 = t->scanner[Acceptor];
				scan3 = t->scanner[Stop];
				score5 = scan5->score(scan5, -exon->end -1);
				score3 = scan3->score(scan3, -exon->start -2);
				break;
			case Exon:
				scan5 = t->scanner[Acceptor];
				scan3 = t->scanner[Donor];
				score5 = scan5->score(scan5, -exon->end -1);
				score3 = scan3->score(scan3, -exon->start +1);
				break;
			default:
				score5 = 0;
				score3 = 0; /* optimizer shush */
				zoeExit("not possible");
		}	
	}
	
	cscore = t->scanner[Coding]->scoref(t->scanner[Coding], exon);
	dscore = zoeScoreDuration(t->hmm->dmap[exon->label], exon->end - exon->start +1);
	xscore = t->exp_score * (exon->end - exon->start +1);
	
	if (error_ok) {
		if (cscore < -1000) cscore = 0;
		if (score5 < -1000) score5 = 0;
		if (score3 < -1000) score3 = 0;
	}
	
	score = score5 + score3 + cscore + dscore - xscore;
	
	if (!padded) {
		exon->start -= PADDING;
		exon->end   -= PADDING;
	}
	
	return score;
}



score_t zoeScoreIntron (zoeTrellis t, zoeFeature intron, int padded) {
	score_t cscore, dscore, xscore;
	
	if (!padded) {
		intron->start += PADDING;
		intron->end   += PADDING;
	}
	
	cscore = t->scanner[Int0]->scoref(t->scanner[Int0], intron);
	dscore = zoeScoreDuration(t->hmm->dmap[Int0], intron->end - intron->start +1);
	xscore = t->exp_score * (intron->end - intron->start +1);
	
	if (!padded) {
		intron->start -= PADDING;
		intron->end   -= PADDING;
	}
	
	return cscore + dscore - xscore;
}


#endif

