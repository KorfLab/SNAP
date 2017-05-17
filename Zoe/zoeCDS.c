/******************************************************************************\
zoeCDS.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_CDS_C
#define ZOE_CDS_C

#include <stdio.h>
#include <string.h>

#include "zoeCDS.h"
#include "zoeFeatureTable.h"

static int MIN_INTRON =     30;
static int MAX_INTRON = 100000;
static int MIN_EXON   =      6;
static int MAX_EXON   =  50000;
static int MIN_GENE   =    150;
static int MAX_GENE   = 500000;
static int MIN_CDS    =    150;
static int MAX_CDS    =  50000;

static const char S5P[5] = {'A', 'C', 'G', 'T', 'N'}; /* positive */
static const char S5N[5] = {'T', 'G', 'C', 'A', 'N'}; /* negative */

/******************************************************************************\
 PRIVATE FUNCTIONS
\******************************************************************************/

static zoeFeatureVec get_sorted_exons (zoeFeatureVec exons, zoeFeatureVec features) {
	int i;
	zoeFeature f;
	
	exons = zoeNewFeatureVec();
	for (i = 0; i < features->size; i++) {
		f = features->elem[i];
		switch(f->label) {
			case Einit: case Eterm: case Exon: case Esngl:
				zoePushFeatureVec(exons, f);
				break;
			default: break;
		}
	}
	
	qsort(exons->elem, exons->size, sizeof(zoeFeature), zoeFeatureCmpPtr);
	
	return exons;
}

static int idiotic_errors (zoeCDS cds) {
	int        i, idx, plus_count = 0, minus_count = 0;
	int        Esngl_count = 0, Einit_count = 0, Eterm_count = 0;
	char       text[256];
	zoeFeature exon;
	
	/* identify strand */
	for (i = 0; i < cds->exons->size; i++) {
		exon = cds->exons->elem[i];
		if (exon->strand == '+') plus_count++;
		if (exon->strand == '-') minus_count++;
	}
	if (plus_count && minus_count) {
		sprintf(text, "gene:mixed_strands");
		zoePushTVec(cds->errors, text);
	} else if (plus_count) {
		cds->strand = '+';
	} else if (minus_count) {
		cds->strand = '-';
	} else {
		sprintf(text, "gene:unknown_strand");
		zoePushTVec(cds->errors, text);
	}
	
	/* change all undefined strands to consensus strand */
	for (i = 0; i < cds->exons->size; i++) {
		exon = cds->exons->elem[i];
		if (exon->strand == UNDEFINED_STRAND) exon->strand = cds->strand;
	}
	
	/* look for multiple occurrances of unique exon types */
	for (i = 0; i < cds->exons->size; i++) {
		exon = cds->exons->elem[i];
		switch (exon->label) {
			case Esngl: Esngl_count++; break;
			case Einit: Einit_count++; break;
			case Eterm: Eterm_count++; break;
			default: break;
		}
	}
	if (Esngl_count > 1) {
		sprintf(text, "gene:multiple_Esngl");
		zoePushTVec(cds->errors, text);
	}
	if (Einit_count > 1) {
		sprintf(text, "gene:multiple_Einit");
		zoePushTVec(cds->errors, text);
	}
	if (Eterm_count > 1) {
		sprintf(text, "gene:multiple_Eterm");
		zoePushTVec(cds->errors, text);
	}
	
	/* make sure Einit is at start and Eterm is at end */
	if (Einit_count) {
		if (cds->strand == '+') {
			if (cds->exons->elem[0]->label != Einit) {
				sprintf(text, "gene:misordered_Einit");
				zoePushTVec(cds->errors, text);
			}
		} else {
			if (cds->exons->elem[cds->exons->size-1]->label != Einit) {
				sprintf(text, "gene:misordered_Einit");
				zoePushTVec(cds->errors, text);
			}
		}
	}
	if (Eterm_count) {
		if (cds->strand == '+') {
			if (cds->exons->elem[cds->exons->size-1]->label != Eterm) {
				sprintf(text, "gene:misordered_Eterm");
				zoePushTVec(cds->errors, text);
			}
		} else {
			if (cds->exons->elem[0]->label != Eterm) {
				sprintf(text, "gene:misordered_Eterm");
				zoePushTVec(cds->errors, text);
			}
		}
	}
	
	/* look for exons out of bounds */
	for (i = 0; i < cds->exons->size; i++) {
		exon = cds->exons->elem[i];
		idx = (cds->strand == '+') ? i +1 : cds->exons->size -i;
		if (exon->start > cds->dna->length -1 || exon->end > cds->dna->length -1) {
			sprintf(text, "exon-%d:out_of_bounds", idx);
			zoePushTVec(cds->errors, text);
		}
	}
	
	/* look for overlapping exons */
	for (i = 1; i < cds->exons->size; i++) {
		exon = cds->exons->elem[i];
		idx = (cds->strand == '+') ? i +1 : cds->exons->size -i;
		if (cds->exons->elem[i-1]->end >= exon->start) {
			sprintf(text, "exon-%d:overlaps_prev_exon", idx);
			zoePushTVec(cds->errors, text);
		} else if (cds->exons->elem[i-1]->end + 1 == exon->start) {
			sprintf(text, "exon-%d:abuts_prev_exon", idx);
			zoePushTVec(cds->errors, text);
		}
	}
		
	if (cds->errors->size) return 1;
	else                   return 0;
}

static zoeDNA transcribe (const char * name, const zoeDNA dna, const zoeFeatureVec exons) {
	zoeFeature   exon = NULL;
	zoeDNA       tx, anti;
	int          i, estart, elength, nt_length = 0;
	char       * seq = NULL;
			
	/* calculate length of transcript */
	for (i = 0; i < exons->size; i++) {
		exon = exons->elem[i];
		nt_length += (exon->end - exon->start + 1);
	}
	
	/* create transcript char[] */
	seq = zoeMalloc(nt_length +1);
	estart = 0;
	for (i = 0; i < exons->size; i++) {
		exon = exons->elem[i];
		elength = exon->end - exon->start + 1;
		memcpy(&seq[estart], &dna->seq[exon->start], elength);
		estart += elength;
	}
	seq[nt_length] = '\0';
	if (strlen(seq) != nt_length) zoeExit("transcribe fatal error");
	
	if (exons->elem[0]->strand == '+') {
		tx = zoeNewDNA(name, seq);
		free(seq);
		return tx;
	} else {
		tx = zoeNewDNA(name, seq);
		anti = zoeAntiDNA(name, tx);
		zoeDeleteDNA(tx);
		free(seq);
		return anti;
	}
}

static int has_internal_stops (zoeProtein pro) {
	int i, count = 0;
	for (i = 0; i < pro->length -1; i++) if (pro->seq[i] == '*') count++;
	return count;
}

/*
static int has_ambiguous_aa (zoeProtein pro) {
	int i, count = 0;
	for (i = 0; i < pro->length; i++) if (pro->seq[i] == 'X') count++;
	return count;
}
*/

struct exon_feature {
	int start;
	int stop;
	int lazy; /* lazy stop - 1 codon away */
	int acceptor;
	int donor;
};

static struct exon_feature get_exon_features (zoeFeature exon, zoeDNA dna) {
	struct exon_feature ef = {0,0,0,0,0};
	char c1, c2, c3;
	
	/* start */
	if (exon->strand == '+') {
		c1 = dna->s5[exon->start];
		c2 = dna->s5[exon->start +1];
		c3 = dna->s5[exon->start +2];
		if (c1 == 0 && c2 == 3 && c3 == 2) ef.start = 1;
	} else if (exon->strand == '-') {
		c1 = dna->s5[exon->end];
		c2 = dna->s5[exon->end -1];
		c3 = dna->s5[exon->end -2];
		if (c1 == 3 && c2 == 0 && c3 == 1) ef.start = 1;
	}
	
	/* stop */
	if (exon->strand == '+') {
		c1 = dna->s5[exon->end -2];
		c2 = dna->s5[exon->end -1];
		c3 = dna->s5[exon->end];
		if (c1 == 3) {
			if      (c2 == 0 && c3 == 0) ef.stop = 1;
			else if (c2 == 0 && c3 == 2) ef.stop = 1;
			else if (c2 == 2 && c3 == 0) ef.stop = 1;
		}
	} else if (exon->strand == '-') {
		c1 = dna->s5[exon->start +2];
		c2 = dna->s5[exon->start +1];
		c3 = dna->s5[exon->start];
		if (c1 == 0) {
			if      (c2 == 3 && c3 == 3) ef.stop = 1;
			else if (c2 == 3 && c3 == 1) ef.stop = 1;
			else if (c2 == 1 && c3 == 3) ef.stop = 1;
		}
	}
	
	/* lazy stop */
	if (exon->strand == '+' && ef.stop == 0 && exon->end +3 < dna->length) {
		c1 = dna->s5[exon->end +1];
		c2 = dna->s5[exon->end +2];
		c3 = dna->s5[exon->end +3];
		if (c1 == 3) {
			if      (c2 == 0 && c3 == 0) ef.lazy = 1;
			else if (c2 == 0 && c3 == 2) ef.lazy = 1;
			else if (c2 == 2 && c3 == 0) ef.lazy = 1;
		}
	} else if (exon->strand == '-' && ef.stop == 0 && exon->start >= 3) {
		c1 = dna->s5[exon->start -1];
		c2 = dna->s5[exon->start -2];
		c3 = dna->s5[exon->start -3];
		if (c1 == 0) {
			if      (c2 == 3 && c3 == 3) ef.lazy = 1;
			else if (c2 == 3 && c3 == 1) ef.lazy = 1;
			else if (c2 == 1 && c3 == 3) ef.lazy = 1;
		}
	}
	
	/* donor */
	if (exon->strand == '+' && exon->end +2 < dna->length) {
		c1 = dna->s5[exon->end+1];
		c2 = dna->s5[exon->end+2];
		if (c1 == 2 && c2 == 3) ef.donor = 1;
	} else if (exon->strand == '-' && exon->start >= 2) {
		c1 = dna->s5[exon->start -1];
		c2 = dna->s5[exon->start -2];
		if (c1 == 1 && c2 == 0) ef.donor = 1;
	}
	
	/* acceptor */
	if (exon->strand == '+' && exon->start >= 2) {
		c1 = dna->s5[exon->start -2];
		c2 = dna->s5[exon->start -1];
		if (c1 == 0 && c2 == 2) ef.acceptor = 1;
	} else if (exon->strand == '-' && exon->end +2 < dna->length) {
		c1 = dna->s5[exon->end +2];
		c2 = dna->s5[exon->end +1];
		if (c1 == 3 && c2 == 1) ef.acceptor = 1;
	}
	
	return ef;
}

struct best_translation {
	strand_t inc5;
	strand_t inc3;
	zoeProtein aa;
};

static struct best_translation get_best_translation (zoeDNA tx) {
	struct best_translation bt;
	zoeProtein pro[3];
	int i, best_score, best_idx, score;
	
	best_score = -1000000;
	best_idx   = -1;
	for (i = 0; i <= 2; i++) {
		pro[i] = zoeTranslateDNA(tx->def, tx, i);
		score = - 3 * has_internal_stops(pro[i]);
		if (pro[i]->seq[pro[i]->length -1] == '*') score++;
		if (pro[i]->seq[0] == 'M') score++;
		if (score > best_score) {
			best_score = score;
			best_idx = i;
		}
	}
	
	bt.inc5 = best_idx;
	bt.inc3 = (tx->length - best_idx) % 3;
	bt.aa   = pro[best_idx];
	
	for (i = 0; i <= 2; i++) if (i != best_idx) zoeDeleteProtein(pro[i]);
	
	return bt;
}

/******************************************************************************\
 PUBLIC FUNCTIONS
\******************************************************************************/

void zoeDeleteCDS (zoeCDS cds) {	
	if (cds == NULL) return;
	
	if (cds->dna) {
		cds->dna = NULL; /* does not delete, just a pointer to parent */
	}
	if (cds->name) {
		zoeFree(cds->name);
		cds->name = NULL;
	}
	if (cds->exons) {
		zoeDeleteFeatureVec(cds->exons);
		cds->exons = NULL;
	}
	if (cds->introns) {
		zoeDeleteFeatureVec(cds->introns);
		cds->introns = NULL;
	}
	if (cds->tx) {
		zoeDeleteDNA(cds->tx);
		cds->tx = NULL;
	}
	if (cds->aa) {
		zoeDeleteProtein(cds->aa);
		cds->aa = NULL;
	}
	if (cds->warnings) {
		zoeDeleteTVec(cds->warnings);
		cds->warnings = NULL;
	}
	if (cds->errors) {
		zoeDeleteTVec(cds->errors);
		cds->errors = NULL;
	}
	zoeFree(cds);
	cds = NULL;
}

zoeCDS zoeNewCDS (const char * name, const zoeDNA dna, const zoeFeatureVec features) {
	zoeCDS cds = zoeMalloc(sizeof(struct zoeCDS));
	zoeFeature exon, einit, eterm, intron;
	char text[256];
	struct exon_feature efstart, efstop;
	struct best_translation bt;	
	int i, nt_length, elength, inc5, inc3, frame, idx, b1, b2, b3, b4, split;
	int exstart, exend, exinc;
	
	/* assigned immediately */
	cds->dna      = dna;
	cds->name     = zoeMalloc(strlen(name) + 1); strcpy(cds->name, name);
	cds->exons    = get_sorted_exons(cds->exons, features);
	cds->source   = zoeCopyFeatureVec(features);
	cds->score    = MIN_SCORE; /* not assigned in zoeNewCDS */
	cds->strand   = '='; /* edited in idiotic_errors */
	
	/* assigned at various points */
	cds->errors   = zoeNewTVec();
	cds->warnings = zoeNewTVec();
	
	/* assigned late */
	cds->introns     = zoeNewFeatureVec();
	cds->start       = 0;
	cds->end         = 0;
	cds->inc5        = 0;
	cds->inc3        = 0;
	cds->tx          = NULL;
	cds->aa          = NULL;
	cds->start_found = 0;
	cds->end_found   = 0;
	cds->OK          = 0;
	
	/* check for no exons - can happen in SNAP that labeled things aren't exons */
	if (cds->exons->size == 0) {
		sprintf(text, "gene:no_exons");
		zoePushTVec(cds->errors, text);
		return cds;
	}
	
	/* check for idiotic errors and abort if necessary */
	if (idiotic_errors(cds)) {
		zoeReportCDS(stderr, cds);
		return cds;
	}
	
	/* start and stop */
	if (cds->strand == '+') {
		einit = cds->exons->elem[0];
		eterm = cds->exons->elem[cds->exons->size -1];
		efstart = get_exon_features(einit, dna);
		efstop = get_exon_features(eterm, dna);
		if (efstop.lazy) eterm->end += 3;
		cds->start = einit->start;
		cds->end   = eterm->end;
		
	} else {
		einit = cds->exons->elem[cds->exons->size -1];
		eterm = cds->exons->elem[0];
		efstart = get_exon_features(einit, dna);		
		efstop = get_exon_features(eterm, dna);
		if (efstop.lazy) eterm->start -= 3;
		cds->start = eterm->start;
		cds->end   = einit->end;
	}
	
	nt_length = 0;
	for (i = 0; i < cds->exons->size; i++) {
		nt_length += cds->exons->elem[i]->end - cds->exons->elem[i]->start +1;
	}
	
	if (efstart.start) cds->start_found = 1;
	if (nt_length % 3 == 0)
		if (efstop.lazy || efstop.stop) cds->end_found = 1;
	
	
	if (!(cds->start_found && cds->end_found)) {
		split = 0;
		if (einit->end - einit->start < 2) {
			sprintf(text, "split-start");
			split = 1;
		}
		if (eterm->end - eterm->start < 2) {
			sprintf(text, "split-stop");
			split = 1;
		}
		if (split == 0) {
			sprintf(text, "cds:incomplete");
		}
		zoePushTVec(cds->warnings, text);
	}
	
	/* re-label generic exons as Einit, Eterm, Esngl only if all exons defined as Exon 
	all_exon = 1;
	for (i = 0; i < cds->exons->size; i++) {
		if (cds->exons->elem[i]->label != Exon) {
			all_exon = 0;
			break;
		}
	}
	*/
	
	/*
	if (all_exon) {
		if (cds->exons->size == 1) {
			exon = cds->exons->elem[0];
			if (cds->start_found && cds->end_found) exon->label = Esngl;
			else if (cds->start_found)              exon->label = Einit;
			else if (cds->end_found)                exon->label = Eterm;
		} else {
			if (cds->start_found) einit->label = Einit;
			if (cds->end_found)   eterm->label = Eterm;
		}
	}
	*/
	
	/* transcribe and translate */
	cds->tx = transcribe(name, dna, cds->exons);
	bt = get_best_translation(cds->tx);
	cds->inc5 = bt.inc5;
	cds->inc3 = bt.inc3;
	cds->aa   = bt.aa;
	
	if (has_internal_stops(cds->aa)) {
		sprintf(text, "cds:internal_stop");
		zoePushTVec(cds->errors, text);
	}
	
	/*
	if (has_ambiguous_aa(cds->aa)) {
		sprintf(text, "cds:ambiguous_aa");
		zoePushTVec(cds->warnings, text);
	}
	*/
	
	/* causing probs
	if (cds->start_found && cds->inc5 != 0) {
		sprintf(text, "cds:start_conflict");
		zoePushTVec(cds->errors, text);
	}
	if (cds->end_found && cds->inc3 != 0) {
		sprintf(text, "cds:stop_conflict");
		zoePushTVec(cds->errors, text);
	}*/
	
	/* create introns */
	intron = zoeNewFeature(Intron, 0, 0, cds->strand, 0, 0, 0, 0, name);
	for (i = 1; i < cds->exons->size; i++) {
		intron->start  = cds->exons->elem[i -1]->end +1;
		intron->end    = cds->exons->elem[i]->start -1;
		zoePushFeatureVec(cds->introns, intron);
	}
	zoeDeleteFeature(intron);
	
	/* short/long genes */
	nt_length = cds->end - cds->start + 1;
	if (nt_length < MIN_GENE) {
		sprintf(text, "gene:short(%d)", nt_length);
		zoePushTVec(cds->warnings, text);
	}
	if (nt_length > MAX_GENE) {
		sprintf(text, "gene:long(%d)", nt_length);
		zoePushTVec(cds->warnings, text);
	}
	
	/* short/long exons */
	for (i = 0; i < cds->exons->size; i++) {
		exon = cds->exons->elem[i];		
		idx = (cds->strand == '+') ? i +1 : cds->exons->size -i;
		nt_length = exon->end - exon->start + 1;
		if (nt_length < MIN_EXON) {
			sprintf(text, "exon-%d:short(%d)", idx, nt_length);
			zoePushTVec(cds->warnings, text);
		} if (nt_length > MAX_EXON) {
			sprintf(text, "exon-%d:long(%d)", idx, nt_length);
			zoePushTVec(cds->warnings, text);
		}
	}
	
	/* short/long introns */
	for (i = 0; i < cds->introns->size; i++) {
		intron = cds->introns->elem[i];
		idx = (cds->strand == '+') ? i +1 : cds->introns->size -i;
		nt_length = intron->end - intron->start + 1;
		if (nt_length < MIN_INTRON) {
			sprintf(text, "intron-%d:short(%d)", idx, nt_length);
			zoePushTVec(cds->warnings, text);
		}
		if (nt_length > MAX_INTRON) {
			sprintf(text, "intron-%d:long(%d)", idx, nt_length);
			zoePushTVec(cds->warnings, text);
		}
	}
	
	/* canonical splicing */
	for (i = 0; i < cds->introns->size; i++) {
		intron = cds->introns->elem[i];
		idx = (cds->strand == '+') ? i +1 : cds->introns->size -i;
		
		b1 = cds->dna->s5[intron->start];
		b2 = cds->dna->s5[intron->start +1];
		b3 = cds->dna->s5[intron->end -1];
		b4 = cds->dna->s5[intron->end];
		
		if (intron->strand == '+') {
			if ( !(b1==2 && b2==3 && b3==0 && b4==2)) {
				sprintf(text, "intron-%d:%c%c..%c%c", idx, S5P[b1], S5P[b2], S5P[b3], S5P[b4]);
				zoePushTVec(cds->warnings, text);
			}
		} else {
			if ( !(b1==1 && b2==3 && b3==0 && b4==1)) {
				sprintf(text, "intron-%d:%c%c..%c%c", idx, S5N[b4], S5N[b3], S5N[b2], S5N[b1]);
				zoePushTVec(cds->warnings, text);
			}
		}
	}
	
	/* label with correct phase and frame */
	if      (cds->inc5 == 1) nt_length = 2;
	else if (cds->inc5 == 2) nt_length = 1;
	else                     nt_length = 0;
		
	elength = 0;
	
	/* dgg; fix for cdsincomplete/complete: phase/inc5 should be 0 for complete */
	if (cds->strand == '-') {
		exstart = cds->exons->size -1;
		exend = -1;
		exinc = -1;
	} else {
		exstart = 0;
		exend = cds->exons->size;
		exinc = 1;
	}
	
	for (i = exstart; i != exend; i += exinc) { /* dgg fixes */
		exon        = cds->exons->elem[i];
		elength     = exon->end - exon->start + 1;
		nt_length  += elength;
		inc3        = nt_length % 3;
		inc5        = (elength - inc3) % 3;
		frame       = (exon->start + inc5) % 3;
		exon->frame = frame;
		exon->inc5  = inc5;
		if (exon->inc5 == -1) exon->inc5 = 2;
		exon->inc3  = inc3;
		if (!zoeVerifyFeature(exon)) {
			zoeWarn("exon does not validate after correcting phase & frame");
			zoeExit("%s", cds->name);
		}
	}
	
	if (cds->errors->size == 0) cds->OK = 1;
	
	return cds;
}

void zoeAntiCDS (zoeCDS cds, int length) {
	int i, start, end;
	
	start = length - cds->end   -1;
	end   = length - cds->start -1;
	
	cds->start  = start;
	cds->end    = end;
	cds->strand = (cds->strand == '+') ? '-' : '+';
	
	for (i = 0; i < cds->exons->size; i++) {
		zoeAntiFeature(cds->exons->elem[i], length);
	}
	
	for (i = 0; i < cds->introns->size; i++) {
		zoeAntiFeature(cds->introns->elem[i], length);
	}
}

void zoeWriteCDS (FILE * stream, const zoeCDS cds) {
	int i;
	
	/*zoeS(stream, ">%s\n", cds->name);*/
	for (i = 0; i < cds->exons->size; i++) {
		zoeWriteFeature(stream, cds->exons->elem[i]);
	}
}

void zoeWriteFullCDS (FILE * stream, const zoeCDS cds) {
	int i;
	
	/*zoeS(stream, ">%s\n", cds->name);*/
	for (i = 0; i < cds->exons->size; i++) {
		zoeWriteFeature(stream, cds->exons->elem[i]);
		if (i < cds->introns->size)
			zoeWriteFeature(stream, cds->introns->elem[i]);
	}
}

void zoeWriteTriteCDS (FILE * stream, const zoeCDS cds) {
	int i;
	
	/*zoeS(stream, ">%s\n", cds->name);*/
	for (i = 0; i < cds->exons->size; i++) {
		zoeWriteTriteFeature(stream, cds->exons->elem[i]);
	}
}

void zoeReportCDS (FILE * stream, const zoeCDS cds) {
	int i;
	
	zoeS(stream, "%s %d %d %d %c", cds->name, cds->start +1, cds->end +1,
		cds->exons->size, cds->strand);
	
	if (cds->errors->size) {
		zoeS(stream, " errors(%d):", cds->errors->size);
		for (i = 0; i < cds->errors->size; i++) {
			zoeS(stream, " %s", cds->errors->elem[i]);
		}
	}
	
	if (cds->warnings->size) {
		zoeS(stream, " warnings(%d):", cds->warnings->size);
		for (i = 0; i < cds->warnings->size; i++) {
			zoeS(stream, " %s", cds->warnings->elem[i]);
		}
	}
	
	zoeS(stream, "\n");
}

int zoeCDScmpptr (const void * v1, const void * v2) {
	return zoeCDScmp( *(zoeCDS *)v1, *(zoeCDS *)v2 );
}

int zoeCDScmp (const zoeCDS f1, const zoeCDS f2) {

	/* coorindate checks: starts then ends */
	if      (f1->start < f2->start) return -1;
	else if (f1->start > f2->start) return  1;
	if      (f1->end < f2->end) return -1;
	else if (f1->end > f2->end) return  1;
	
	/* strand check */
	if      (f1->strand < f2->strand) return -1;
	else if (f1->strand > f2->strand) return  1;

	/* transcript length check */
	if (f1->tx && f2->tx) {
		if      (f1->tx->length > f2->tx->length) return  1;
		else if (f1->tx->length < f2->tx->length) return -1;
	}
	
	/* transcript sequence check */
	if (f1->tx && f2->tx) {
		if (strcmp(f1->tx->seq, f2->tx->seq) == 0) return 0;
		else                                       return 1;
	}
	
	/* I suppose I could check the individual features ... */
	return 0;

}

int zoeCDSsOverlap (const zoeCDS a, const zoeCDS b) {
	if      (a->start >= b->start && a->start <= b->end) return 1;
	else if (b->start >= a->start && b->start <= a->end) return 1;
	else                                                 return 0;
}

int zoeCDSsShareSequence (const zoeCDS a, const zoeCDS b) {
	int i, j;
	
	if (!zoeCDSsOverlap(a, b)) return 0;
	 
	for (i = 0; i < a->exons->size; i++) {
		for (j = 0; j < b->exons->size; j++) {
			if (zoeFeaturesOverlap(a->exons->elem[i], b->exons->elem[j])) {
				return 1;
			}
		}
	}
	
	return 0;
}

void zoeSetMinIntron (int i) {MIN_INTRON = i;}
void zoeSetMaxIntron (int i) {MAX_INTRON = i;}
void zoeSetMinExon (int i)   {MIN_EXON   = i;}
void zoeSetMaxExon (int i)   {MAX_EXON   = i;}
void zoeSetMinGene (int i)   {MIN_GENE   = i;}
void zoeSetMaxGene (int i)   {MAX_GENE   = i;}
void zoeSetMinCDS (int i)    {MIN_CDS    = i;}
void zoeSetMaxCDS (int i)    {MAX_CDS    = i;}

#endif
