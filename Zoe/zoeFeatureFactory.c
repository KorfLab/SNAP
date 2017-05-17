/******************************************************************************\
 zoeFeatureFactory.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FEATUREFACTORY_C
#define ZOE_FEATUREFACTORY_C

#include "zoeFeatureFactory.h"

static const int PADDING = 48; /* as in trellis, should make it extern */

/******************************************************************************\
 PRIVATE FUNCTIONS
\******************************************************************************/

static zoeFeatureVec zoeMakeFeatures (const zoeFeatureFactory fac, coor_t pos) {
	score_t       score;
	zoeFeature    f;
	zoeFeatureVec vec;

	score = fac->scanner->score(fac->scanner, pos);
	if (score == MIN_SCORE) return NULL;
	
	f = zoeNewFeature(fac->type, pos, pos, '+', score, 0, 0, 0, NULL/*, NULL*/);
	vec = zoeNewFeatureVec();
	zoePushFeatureVec(vec, f);
	zoeDeleteFeature(f);
	return vec;
}

static zoeFeatureVec zoeMakeRepeats (const zoeFeatureFactory fac, coor_t pos) {
	coor_t        i, start, end;
	zoeFeature    f;    /* = zoeNewFeature(None, 0, 0, '+', 0, 0, 0, 0, NULL);*/
	zoeFeatureVec fvec; /* = zoeNewFeatureVec(); */

	/* does a repeat end here? */
	if (fac->dna->s5[pos] == 4) {
		if (pos == fac->dna->length -1) {
			end = fac->dna->length -1;
		} else if (fac->dna->s5[pos +1] != 4) {
			end = pos;
		} else {
			return NULL;
		}
	} else {
		return NULL;
	}
		
	/* find repeat start */
	i = pos;
	while (i > 0) {
		if (fac->dna->s5[i] == 4) i--;
		else                     break;
	}
	start = i;
	
	/* min length filtering */
	if (end - start + 1 < fac->length) return NULL;
	
	/* make a repeat */
	fvec = zoeNewFeatureVec();
	f = zoeNewFeature(Repeat, start, end, '=', 0, 0, 0, 0, NULL/*, NULL*/);
	f->score = fac->scanner->scoref(fac->scanner, f);
	zoePushFeatureVec(fvec, f);	
	zoeDeleteFeature(f);
	return fvec;
}

static zoeFeatureVec zoeMakeExons (const zoeFeatureFactory fac, coor_t pos) {
	int           i, frame = -1, inc5, inc3, begin, end, length;
	zoeFeature    exon = zoeNewFeature(None, 0, 0, '+', 0, 0, 0, 0, NULL/*, NULL*/);
	zoeFeature    e = NULL;
	zoeFeatureVec vec = zoeNewFeatureVec();
	/*int j;*/
	
	/* invariant */
	exon->group  = NULL;
	exon->strand = '+';
	
	/* Esngl */
	if (fac->stop[pos] != MIN_SCORE) {
		frame = pos % 3; /* Esngl is in the same frame as the stop codon */
		end   = pos -3;
		begin = fac->fstop[end];
		for (i = begin; i < end; i += 3) {
			if (fac->start[i] == MIN_SCORE) continue;
			exon->label = Esngl;
			exon->start = i;
			exon->end   = pos -1;
			exon->score = fac->start[i] + fac->stop[pos];
			exon->inc5   = 0;
			exon->inc3   = 0;
			exon->frame = frame;
			exon->group = NULL;
			zoePushFeatureVec(vec, exon);
		}
	}
	
	/* Eterm */
	if (fac->stop[pos] != MIN_SCORE) {
		frame = pos % 3; /* Eterm frame is the same as stop */
		end   = pos -3;
		begin = fac->fstop[end];
		for (i = begin; i < end; i++) {
			if (fac->acc[i] == MIN_SCORE) continue;
			exon->label = Eterm;
			exon->start = i + 1;
			exon->end   = pos -1;
			exon->score = fac->acc[i] + fac->stop[pos];
			exon->inc5   = (exon->end - exon->start + 1) % 3;
			exon->inc3   = 0;
			exon->frame = frame % 3;
			zoePushFeatureVec(vec, exon);
		}
	}
	
	/* Einit */
	if (fac->don[pos] != MIN_SCORE) {		
		for (inc3 = 0; inc3 < 3; inc3++) {
			frame = (pos - inc3) % 3;
			end   = pos - 3;
			begin = fac->fstop[end -inc3];
			
			/* for (i = begin; i < end; i += 3) { - old code */
			/* for (i = begin; i <= end; i += 3) { - new code */
			
			for (i = begin; i < end; i += 3) {
				if (fac->start[i] == MIN_SCORE) continue;
				exon->label = Einit;
				exon->start = i;
				exon->end   = pos - 1;
				exon->score = fac->start[i] + fac->don[pos];
				exon->inc5   = 0;
				exon->inc3   = inc3;
				exon->frame = frame;
				zoePushFeatureVec(vec, exon);
			}
		}
	}
		
	/* Exon (internal) */
	if (fac->don[pos] != MIN_SCORE) {
		for (inc3 = 0; inc3 < 3; inc3++) {
			frame = (pos - inc3) % 3;
			end   = pos -3;
			begin = fac->fstop[end -inc3];
			for (i = begin; i < end; i ++) {
				if (fac->acc[i] == MIN_SCORE) continue;
				length = pos -1 - i;
				inc5 = (length - inc3) % 3;
				exon->label  = Exon;
				exon->start  = i + 1;
				exon->end    = pos - 1;
				exon->score  = fac->acc[i] + fac->don[pos];
				exon->inc5    = inc5;
				exon->inc3    = inc3;
				exon->frame  = frame;
				zoePushFeatureVec(vec, exon);
			}
		}
	}
			
	/* add coding score */
	for (i = 0; i < vec->size; i++) {
		e = vec->elem[i];
		
		/* convert exon frame to cds[frame] -  I know, it looks weird, trust me... */
		switch (e->frame) {
			case 0: frame = 2; break;
			case 1: frame = 1; break;
			case 2: frame = 0; break;
			default: zoeExit("zoeMakeExons impossible");
		}
		
		/*
		HARD-CODED FEATURE: Since acceptor and donor overlap each exon
		by 3 bp (by my convention), these bases must not be scored by
		the coding model
		
		--------[===EXON===]--------
		 <----------  Acceptor overlaps by 3
		                ---------> Donor overlaps by 3
		*/
		
		length = e->end - e->start + 1;
		if (length > 12) {
			e->score += fac->cds[frame][e->end -3]
				      - fac->cds[frame][e->start +3]; /* changed +9 to +3 */
		}
	}
		
	zoeDeleteFeature(exon);
	
	if (vec->size) {
		return vec;
	} else {
		zoeDeleteFeatureVec(vec);
		return NULL;
	}
}

static zoeFeatureVec zoeMakeOpenReadingFrames (const zoeFeatureFactory fac, coor_t pos) {
	zoeFeatureVec vec;
	zoeFeature    orf;
	char          key[16];
	
	sprintf(key, "%d", pos);
	orf = zoeGetHash(fac->hash, key);
	if (orf == NULL) return NULL;
	
	vec = zoeNewFeatureVec();
	zoePushFeatureVec(vec, orf);
	return vec;
}

/******************************************************************************\
 PUBLIC FUNCTIONS
\******************************************************************************/

void zoeDeleteFeatureFactory(zoeFeatureFactory f) {
	int i;
	
	if (f == NULL) return;
	
	if (f->orfs) {
		zoeDeleteFeatureVec(f->orfs);
		f->orfs = NULL;
	}
	
	if (f->hash) {
		zoeDeleteHash(f->hash);
		f->hash = NULL;
	}
	
	for (i = 0; i < 3; i++) {
		if (f->cds[i]) {
			zoeFree(f->cds[i]);
			f->cds[i] = NULL;
		}
	}
	if (f->start) {zoeFree(f->start); f->start = NULL;}
	if (f->stop)  {zoeFree(f->stop);  f->stop  = NULL;}
	if (f->acc)   {zoeFree(f->acc);   f->acc   = NULL;}
	if (f->don)   {zoeFree(f->don);   f->don   = NULL;}
	if (f->fstop) {zoeFree(f->fstop); f->fstop = NULL;}
	zoeFree(f);
	f = NULL;
}

zoeFeatureFactory zoeNewEFactory (
	zoeScanner cds_scan,
	zoeScanner accpt_scan,
	zoeScanner donor_scan,
	zoeScanner start_scan,
	zoeScanner stop_scan)
{

	score_t         * cds[3];
	score_t         * acc;
	score_t         * don;
	score_t         * start;
	score_t         * stop;
	int               f0, f1, f2;
	int             * fstop;	
	coor_t            i;
	zoeDNA            dna = cds_scan->dna;
	zoeFeatureFactory factory = zoeMalloc(sizeof(struct zoeFeatureFactory));
		
	/* cds frame-specific scanners */
	zoeScanner cscan[3];
	cscan[0] = cds_scan->subscanner[0];
	cscan[1] = cds_scan->subscanner[1];
	cscan[2] = cds_scan->subscanner[2];

	/* allocate storage for factory attributes */
	cds[0] = zoeMalloc(dna->length * sizeof(score_t));
	cds[1] = zoeMalloc(dna->length * sizeof(score_t));
	cds[2] = zoeMalloc(dna->length * sizeof(score_t));
	acc    = zoeMalloc(dna->length * sizeof(score_t));	
	don    = zoeMalloc(dna->length * sizeof(score_t));
	start  = zoeMalloc(dna->length * sizeof(score_t));
	stop   = zoeMalloc(dna->length * sizeof(score_t));
	fstop  = zoeMalloc(dna->length * sizeof(int));
	
	/* compute feature positions -------------------------------------------- */
	for (i = 0; i < dna->length; i++) {
		acc[i]   = accpt_scan->score(accpt_scan, i);
		don[i]   = donor_scan->score(donor_scan, i);
		start[i] = start_scan->score(start_scan, i);
		stop[i]  = stop_scan->score(stop_scan, i);
	}
	
	/* compute CDS scores in 3 frames --------------------------------------- */
	cds[0][0] = 0;
	cds[1][0] = 0;
	cds[2][0] = 0;
	for (i = 1; i < dna->length; i++) {
		f0 = (i - 1) % 3;
		f1 = (i + 0) % 3;
		f2 = (i + 1) % 3;
		cds[0][i] = cscan[f0]->score(cscan[f0], i);
		cds[1][i] = cscan[f1]->score(cscan[f1], i);
		cds[2][i] = cscan[f2]->score(cscan[f2], i);
		
		if (cds[0][i] == MIN_SCORE) cds[0][i]  = 0;
		else                               cds[0][i] += cds[0][i-1];
		if (cds[1][i] == MIN_SCORE) cds[1][i]  = 0;
		else                               cds[1][i] += cds[1][i-1];
		if (cds[2][i] == MIN_SCORE) cds[2][i]  = 0;
		else                               cds[2][i] += cds[2][i-1];
	}
	
	/* positions of stop codons in each frame  ------------------------------ */	
	fstop[0] = 0;
	fstop[1] = 1;
	fstop[2] = 2;
	
	for (i = 3; i < dna->length; i++) {
		if (stop[i] != MIN_SCORE) {
			fstop[i] = i;
		} else {
			fstop[i] = fstop[i-3];
		}
	}

	/* set factory attributes ----------------------------------------------- */
	factory->create  = zoeMakeExons;
	factory->type    = Exon;
	factory->dna     = dna;
	factory->cds[0]  = cds[0];
	factory->cds[1]  = cds[1];
	factory->cds[2]  = cds[2];
	factory->acc     = acc;
	factory->don     = don;
	factory->start   = start;
	factory->stop    = stop;
	factory->fstop   = fstop;
	factory->offset  = cscan[0]->model->focus;
		
	/* this stuff isn't used by an EFactory */
	factory->length  = 0;
	factory->score   = 0;
	factory->scanner = NULL;
	factory->orfs    = NULL;
	factory->hash    = NULL;
	
	return factory;
}

zoeFeatureFactory zoeNewXFactory  (
	zoeScanner cds_scan,
	zoeScanner accpt_scan,
	zoeScanner donor_scan,
	zoeScanner start_scan,
	zoeScanner stop_scan,
	coor_t     min_length,
	score_t    min_score) {

	zoeFeatureFactory factory = zoeMalloc(sizeof(struct zoeFeatureFactory));
	zoeFeatureFactory efac;
	zoeFeatureVec     exons;
	zoeFeature        exon, max_exon;
	score_t           max_score;
	zoeScanner        cs, as, ds, ms, ss;
	int               i, j, length;
	char              key[16];
	
	factory->create  = zoeMakeOpenReadingFrames;
	factory->type    = ORF;
	factory->length  = min_length;
	factory->score   = min_score;
	
	/* not used by XFactory */
	factory->scanner   = NULL;
	factory->strand    = UNDEFINED_STRAND;
	factory->cds[0]    = NULL;
	factory->cds[1]    = NULL;
	factory->cds[2]    = NULL;
	factory->acc       = NULL;
	factory->don       = NULL;
	factory->start     = NULL;
	factory->stop      = NULL;
	factory->fstop     = NULL;
	
	cs = zoeNewScanner(cds_scan->anti, cds_scan->dna, cds_scan->model);
	as = zoeNewScanner(cds_scan->anti, cds_scan->dna, accpt_scan->model);
	ds = zoeNewScanner(cds_scan->anti, cds_scan->dna, donor_scan->model);
	ms = zoeNewScanner(cds_scan->anti, cds_scan->dna, start_scan->model);
	ss = zoeNewScanner(cds_scan->anti, cds_scan->dna, stop_scan->model);
	
	efac = zoeNewEFactory(cs, as, ds, ms, ss);
		
	/* XFactory precomputes all ORFs */
	factory->orfs = zoeNewFeatureVec();
	for (i = PADDING; i < cds_scan->dna->length -PADDING; i++) {
		exons = efac->create(efac, i);
		if (exons == NULL) continue;
		
		/* find maximum-scoring exon at each position */
		max_score = MIN_SCORE;
		max_exon = NULL;
		for (j = 0; j < exons->size; j++) {
			exon = exons->elem[j];
			if (exon->score < min_score) continue;
			length = exon->end - exon->start + 1;
			if (length < min_length) continue;
			if (exon->score > max_score) {
				max_score = exon->score;
				max_exon  = exon;
			}
		}
		
		if (max_score != MIN_SCORE) {
			max_exon->label = ORF;
			zoeAntiFeature(max_exon, cds_scan->dna->length);
			zoePushFeatureVec(factory->orfs, max_exon);
		}
		
		zoeDeleteFeatureVec(exons);
	}
	
	/* create a lookup */
	factory->hash = zoeNewHash();
	for (i = 0; i < factory->orfs->size; i++) {
		exon = factory->orfs->elem[i];		
		sprintf(key, "%d", exon->end);
		zoeSetHash(factory->hash, key, exon);
	}
	
	zoeDeleteScanner(cs);
	zoeDeleteScanner(as);
	zoeDeleteScanner(ds);
	zoeDeleteScanner(ms);
	zoeDeleteScanner(ss);
	zoeDeleteFeatureFactory(efac);
		
	return factory;
}

zoeFeatureFactory zoeNewSFactory(
	zoeScanner scanner,
	zoeLabel   type)
{
	zoeFeatureFactory factory = zoeMalloc(sizeof(struct zoeFeatureFactory));
	factory->create  = zoeMakeFeatures;
	factory->type    = type;
	factory->dna     = scanner->dna;
	factory->scanner = scanner;
	
	/* not used by SFactory */
	factory->length    = 0;
	factory->score     = 0;
	factory->orfs      = NULL;
	factory->hash      = NULL;
	factory->cds[0]    = NULL;
	factory->cds[1]    = NULL;
	factory->cds[2]    = NULL;
	factory->acc       = NULL;
	factory->don       = NULL;
	factory->start     = NULL;
	factory->stop      = NULL;
	factory->fstop     = NULL;
	
	return factory;
}

zoeFeatureFactory zoeNewRFactory(
	zoeScanner scanner,
	coor_t  length)
{
	zoeFeatureFactory factory = zoeMalloc(sizeof(struct zoeFeatureFactory));
	factory->create  = zoeMakeRepeats;
	factory->type    = Repeat;
	factory->dna     = scanner->dna;
	factory->length  = length;
	factory->scanner = scanner;

	/* not used by RFactory */
	factory->score     = 0;
	factory->orfs      = NULL;
	factory->hash      = NULL;
	factory->cds[0]    = NULL;
	factory->cds[1]    = NULL;
	factory->cds[2]    = NULL;
	factory->acc       = NULL;
	factory->don       = NULL;
	factory->start     = NULL;
	factory->stop      = NULL;
	factory->fstop     = NULL;
		
	return factory;
}

#endif
