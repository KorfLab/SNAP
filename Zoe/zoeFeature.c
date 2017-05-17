/******************************************************************************\
zoeFeature.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FEATURE_C
#define ZOE_FEATURE_C

#include "zoeFeature.h"

void zoeDeleteFeature (zoeFeature f) {	
	if (f == NULL) return;
	if (f->group) {
		zoeFree(f->group);
		f->group = NULL;
	}
	/*
	if (f->subfeatures) {
		zoeDeleteFeatureVec(f->subfeatures);
		f->subfeatures = NULL;
	}
	*/
	zoeFree(f);
	f = NULL;
}

zoeFeature zoeNewFeature (
		zoeLabel     label,
		coor_t       start,
		coor_t       end,
		strand_t     strand,
		score_t      score,
		frame_t      inc5,
		frame_t      inc3,
		frame_t      frame,
		const char * group)
{
	zoeFeature f = zoeMalloc(sizeof(struct zoeFeature));
	
	f->label  = label;
	f->start  = start;
	f->end    = end;
	f->strand = strand;
	f->score  = score;
	f->inc5   = inc5;
	f->inc3   = inc3;
	f->frame  = frame;
	
	if (group) {
		f->group = zoeMalloc(strlen(group) +1);
		strcpy(f->group, group);
	} else {
		f->group = NULL;
	}
	
	if (!zoeVerifyFeature(f)) {
		zoeWriteFeature(stderr, f);
		zoeExit("zoeNewFeature illegal feature");
		return NULL;
	}
		
	return f;
}

zoeFeature zoeNewTriteFeature (zoeLabel label, coor_t start, coor_t end, const char * group) {	
	if (start <= end) {
		return zoeNewFeature(label, start, end, '+', MIN_SCORE, UNDEFINED_FRAME,
			UNDEFINED_FRAME, UNDEFINED_FRAME, group/*, NULL*/);
	} else {
		return zoeNewFeature(label, end, start, '-', MIN_SCORE, UNDEFINED_FRAME,
			UNDEFINED_FRAME, UNDEFINED_FRAME, group/*, NULL*/);
	}
	
}

zoeFeature zoeReadFeature (FILE * stream) {
	char      line[256], Label[64], Start[64], End[64], Strand[64], Score[64],
	          Loh[64], Roh[64], Frame[64], Group[64];
	char    * group = NULL;
	int       short_form;
	coor_t    tmp;
	zoeLabel  label;
	coor_t    start, end;
	strand_t  strand;
	score_t   score;
	frame_t   inc5, inc3, frame;

	if (fgets(line, sizeof(line), stream) == NULL) return NULL;
	
	if (sscanf(line, "%s %s %s %s %s %s %s %s %s", /* group */
			Label, Start, End, Strand, Score, Loh, Roh, Frame, Group) == 9) {
		short_form = 0;
		group      = Group;
	} else if (sscanf(line, "%s %s %s %s %s %s %s %s", /* no group */
			Label, Start, End, Strand, Score, Loh, Roh, Frame) == 8) {
		short_form = 0;
		group      = NULL;
	} else if (sscanf(line, "%s %s %s %s", /* short form w/ group */
			Label, Start, End, Group) == 4) {
		short_form = 1;
		group      = Group;
	} else if (sscanf(line, "%s %s %s", /* short form w/o group */
			Label, Start, End) == 3) {
		short_form = 1;
		group      = NULL;
	} else {
		zoeWarn("zoeReadFeature format not valid");
		return NULL;
	}
		
	label = zoeText2Label(Label);
	if (label == None) {
		zoeWarn("zoeReadFeature: the following line gave me trouble");
		zoeE("%s", line);
	}
	
	if (short_form) {
		start = zoeText2Coor(Start);
		end   = zoeText2Coor(End);
		if (start < end) {
			strand = '+';
		} else if (start > end) {
			strand = '-';
			tmp    = start;
			start  = end;
			end    = tmp;
		} else {
			strand = '=';
		}
		score = MIN_SCORE;
		inc5  = UNDEFINED_FRAME;
		inc3  = UNDEFINED_FRAME;
		frame = UNDEFINED_FRAME;
	} else {
		start  = zoeText2Coor(Start);
		end    = zoeText2Coor(End);
		strand = zoeText2Strand(Strand);
		score  = zoeText2Score(Score);
		inc5   = zoeText2Frame(Loh);
		inc3   = zoeText2Frame(Roh);
		frame  = zoeText2Frame(Frame);
	}
	
	/* subtract 1 from start and end, zoe uses 0-based coordinates internally */
	start--;
	end--;
	
	return zoeNewFeature(label, start, end, strand, score, inc5, inc3, frame, group);
}

zoeFeature zoeReadGFF (FILE * stream) {
	char       line[1024], Source[64], ID[64], Label[64],
	           Start[16], End[16], Strand[8], Score[32], Frame[8], Group[512];
	char     * group;
	zoeLabel   label;
	coor_t     start, end;
	strand_t   strand;
	score_t    score;
	frame_t    inc5, inc3, frame;

	if (fgets(line, sizeof(line), stream) == NULL) return NULL;
	
	if (sscanf(line, "%s %s %s %s %s %s %s %s %s", /* group */
			Source, ID, Label, Start, End, Score, Strand, Frame, Group) == 9) {
		group = Group;
	} else if (sscanf(line, "%s %s %s %s %s %s %s %s", /* no group */
			Source, ID, Label, Start, End, Score, Strand, Frame) == 8) {
		group = NULL;
	} else {
		zoeWarn("zoeReadGFF format not valid");
		return NULL;
	}
	
	label = zoeText2Label(Label);
	start  = zoeText2Coor(Start);
	end    = zoeText2Coor(End);
	strand = zoeText2Strand(Strand);
	score  = zoeText2Score(Score);
	inc5    = zoeText2Frame(Frame); /* inc5 == gff frame-phase */
	inc3    = UNDEFINED_FRAME;
	frame  = UNDEFINED_FRAME;
	
	/* subtract 1 from start and end, zoe uses 0-based coordinates internally */
	start--;
	end--;
	
	return zoeNewFeature(label, start, end, strand, score, inc5, inc3, frame, group/*, NULL*/);
}

void zoeWriteFeature (FILE * stream, const zoeFeature f) {
	char label[16], start[16], end[16], strand[8], score[32],
	     left[8], right[8], frame[8];

	zoeLabel2Text(f->label, label);
	zoeCoor2Text(f->start +1, start);  /* note: must add 1 to start and end because */
	zoeCoor2Text(f->end +1, end);      /* people use 1-based coordinates */
	zoeStrand2Text(f->strand, strand);
	zoeScore2Text(f->score, score);
	zoeFrame2Text(f->inc5, left);
	zoeFrame2Text(f->inc3, right);
	zoeFrame2Text(f->frame, frame);
		
	zoeS(stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			label, start, end, strand, score, left, right, frame);
	
	if (f->group == NULL) zoeS(stream, "\n");
	else                  zoeS(stream, "\t%s\n", f->group);
}

void zoeWriteDebugFeature (FILE * stream, const zoeFeature f) {
	char label[16], start[16], end[16], strand[8], score[32],
	     left[8], right[8], frame[8];

	zoeLabel2Text(f->label, label);
	zoeCoor2Text(f->start +1 -48, start); 
	zoeCoor2Text(f->end +1 -48, end);      
	zoeStrand2Text(f->strand, strand);
	zoeScore2Text(f->score, score);
	zoeFrame2Text(f->inc5, left);
	zoeFrame2Text(f->inc3, right);
	zoeFrame2Text(f->frame, frame);
		
	zoeS(stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			label, start, end, strand, score, left, right, frame);
	
	if (f->group == NULL) zoeS(stream, "\n");
	else                  zoeS(stream, "\t%s\n", f->group);
}


void zoeWriteTriteFeature (FILE * stream, const zoeFeature f) {
	char label[16], start[16], end[16];

	zoeLabel2Text(f->label, label);
	zoeCoor2Text(f->start +1, start);  /* note: must add 1 to start and end because */
	zoeCoor2Text(f->end +1, end);      /* people use 1-based coordinates */
	
	if (f->strand == '-') zoeS(stream, "%s\t%s\t%s", label, end, start);
	else                  zoeS(stream, "%s\t%s\t%s", label, start, end);
	
	if (f->group == NULL) zoeS(stream, "\n");
	else                  zoeS(stream, "\t%s\n", f->group);
}

void zoeWriteGFF (FILE *stream, const zoeFeature f, const char * id, const char * source) {
	char label[16], start[16], end[16], strand[8], score[32];

	zoeLabel2Text(f->label, label);
	zoeCoor2Text(f->start +1, start);  /* note: must add 1 to start and end because */
	zoeCoor2Text(f->end +1, end);      /* people use 1-based coordinates */
	zoeStrand2Text(f->strand, strand);
	zoeScore2Text(f->score, score);
	
	zoeS(stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.", /* frame undefined */
		source, id, label, start, end, score, strand);
	if (f->group == NULL) zoeS(stream, "\n");
	else                  zoeS(stream, "\t%s\n", f->group);
}

int zoeVerifyFeature (const zoeFeature f) {
	int retval = 1;
	
	if (f->start < 0 || f->end < 0) {
		zoeWarn("zoeVerifyFeature: coordinate not positive", f->start);
		retval = 0;
	}
	
	if (f->start > f->end) {
		zoeWarn("zoeVerifyFeature: start > end");
		retval = 0;
	}
	
	if (f->start == UNDEFINED_COOR || f->end == UNDEFINED_COOR) {
		zoeWarn("zoeVerifyFeature: undefined coordinates");
		retval = 0;
	}
	
	if (f->strand != '-' && f->strand != '+' && f->strand != '=') {
		zoeWarn("zoeVerifyFeature: strand must be {+, -, =} (%c)", f->strand);
		retval = 0;
	}
	
	if (f->inc5 != UNDEFINED_FRAME && (f->inc5 < 0 || f->inc5 > 2)) {
		zoeWarn("zoeVerifyFeature: inc5 out of bounds (%d)", f->inc5);
		retval = 0;
	}
	
	if (f->inc3 != UNDEFINED_FRAME && (f->inc3 < 0 || f->inc3 > 2)) {
		zoeWarn("zoeVerifyFeature: inc3 out of bounds (%d)", f->inc3);
		retval = 0;
	}
	
	if (f->frame != UNDEFINED_FRAME && (f->frame < 0 || f->frame > 2)) {
		zoeWarn("zoeVerifyFeature: frame out of bounds (%d)", f->frame);
		retval = 0;
	}
	
	/* passed the tests */
	return retval;
}

zoeFeature zoeCopyFeature (const zoeFeature f) {
	return zoeNewFeature(f->label, f->start, f->end, f->strand, f->score,
		f->inc5, f->inc3, f->frame, f->group);
}

void zoeAntiFeature (zoeFeature f, int length) {
	coor_t   start, end;
	
	start = length - f->end   -1;
	end   = length - f->start -1;
	f->start = start;
	f->end   = end;
	if      (f->strand == '+') f->strand = '-';
	else if (f->strand == '-') f->strand = '+';
	
}

int zoeFeatureCmp (const zoeFeature f1, const zoeFeature f2) {

	if      (f1->start < f2->start) return -1;
	else if (f1->start > f2->start) return 1;
	else {
		if      (f1->end < f2->end) return -1;
		else if (f1->end > f2->end) return 1;
		else {
			if      (f1->strand < f2->strand) return -1;
			else if (f1->strand > f2->strand) return 1;
			else return f1->label - f2->label;
		}
	}
}

int zoeFeatureCmpPtr (const void * v1, const void * v2) {
	return zoeFeatureCmp( *(zoeFeature *)v1, *(zoeFeature *)v2 );
}

int zoeFeaturesOverlap (const zoeFeature f1, const zoeFeature f2) {
	if      (f1->start >= f2->start && f1->start <= f2->end) return 1;
	else if (f2->start >= f1->start && f2->start <= f1->end) return 1;
	else                                                     return 0;
}


zoeLabel zoeText2Label (const char * string) {
	if      (strcmp(string, "None") == 0)       return None;
	
	else if (strcmp(string, "Inter") == 0)      return Inter;
	else if (strcmp(string, "Int0") == 0)       return Int0;
	else if (strcmp(string, "Int1") == 0)       return Int1;
	else if (strcmp(string, "Int1T") == 0)      return Int1T;
	else if (strcmp(string, "Int2") == 0)       return Int2;
	else if (strcmp(string, "Int2TA") == 0)     return Int2TA;
	else if (strcmp(string, "Int2TG") == 0)     return Int2TG;
	else if (strcmp(string, "Intron") == 0)     return Intron;
	else if (strcmp(string, "UTR5") == 0)       return UTR5;
	else if (strcmp(string, "UTR3") == 0)       return UTR3;
	else if (strcmp(string, "Esngl") == 0)      return Esngl;
	else if (strcmp(string, "Einit") == 0)      return Einit;
	else if (strcmp(string, "Eterm") == 0)      return Eterm;
	else if (strcmp(string, "Exon") == 0)       return Exon;
	else if (strcmp(string, "Coding") == 0)     return Coding;
	else if (strcmp(string, "Gene") == 0)       return Gene;
	else if (strcmp(string, "Acceptor") == 0)   return Acceptor;
	else if (strcmp(string, "Donor") == 0)      return Donor;
	else if (strcmp(string, "Start") == 0)      return Start;
	else if (strcmp(string, "Stop") == 0)       return Stop;
	else if (strcmp(string, "Repeat") == 0)     return Repeat;
	else if (strcmp(string, "CNS") == 0)        return CNS;
	else if (strcmp(string, "ORF") == 0)        return ORF;
	else if (strcmp(string, "PolyA") == 0)      return PolyA;
	else if (strcmp(string, "Prom") == 0)       return Prom;
	else if (strcmp(string, "BPS") == 0)        return BPS;
	else if (strcmp(string, "TSS") == 0)        return TSS;
	else if (strcmp(string, "Misc") == 0)       return Misc;
	else if (strcmp(string, "HSP_NN") == 0)     return HSP_NN;
	else if (strcmp(string, "HSP_NA") == 0)     return HSP_NA;
	else if (strcmp(string, "HSP_AN") == 0)     return HSP_AN;
	else if (strcmp(string, "HSP_AA") == 0)     return HSP_AA;
	else {
		zoeWarn("zoeText2Label illegal label (%s)", string);
		return None;
	}
}

void zoeLabel2Text (zoeLabel label, char * string) {
	switch (label) {
		case None:      strcpy(string, "None");       break;
		case Inter:     strcpy(string, "Inter");      break;
		case Int0:      strcpy(string, "Int0");       break;
		case Int1:      strcpy(string, "Int1");       break;
		case Int1T:     strcpy(string, "Int1T");      break;
		case Int2:      strcpy(string, "Int2");       break;
		case Int2TA:    strcpy(string, "Int2TA");     break;
		case Int2TG:    strcpy(string, "Int2TG");     break;
		case Intron:    strcpy(string, "Intron");     break;
		case UTR5:      strcpy(string, "UTR5");       break;
		case UTR3:      strcpy(string, "UTR3");       break;
		case Esngl:     strcpy(string, "Esngl");      break;
		case Einit:     strcpy(string, "Einit");      break;
		case Eterm:     strcpy(string, "Eterm");      break;
		case Exon:      strcpy(string, "Exon");       break;
		case Coding:    strcpy(string, "Coding");     break;
		case Gene:      strcpy(string, "Gene");       break;
		case Acceptor:  strcpy(string, "Acceptor");   break;
		case Donor:     strcpy(string, "Donor");      break;
		case Start:     strcpy(string, "Start");      break;
		case Stop:      strcpy(string, "Stop");       break;
		case Repeat:    strcpy(string, "Repeat");     break;
		case CNS:       strcpy(string, "CNS");        break;
		case ORF:       strcpy(string, "ORF");        break;
		case PolyA:     strcpy(string, "PolyA");      break;
		case Prom:      strcpy(string, "Prom");       break;
		case BPS:       strcpy(string, "BPS");        break;
		case TSS:       strcpy(string, "TSS");        break;
		case Misc:      strcpy(string, "Misc");       break;
		case HSP_NN:    strcpy(string, "HSP_NN");     break;
		case HSP_NA:    strcpy(string, "HSP_NA");     break;
		case HSP_AN:    strcpy(string, "HSP_AN");     break;
		case HSP_AA:    strcpy(string, "HSP_AA");     break;
		default:
			zoeWarn("zoeLabel2Text illegal feature name, assigning Unknown");
			(void)strcpy(string, "Unknown");
	}
}

void zoeWriteLabel (FILE * stream, zoeLabel label) {
	char string[32];
	zoeLabel2Text(label, string);
	zoeS(stream, string);
}

/* zoeFeatureVec stuff */

void zoeDeleteFeatureVec (zoeFeatureVec vec) {
	int i;
		
	if (vec->elem) {
		for (i = 0; i < vec->size; i++) {
			zoeDeleteFeature(vec->elem[i]);
		}
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeFeatureVec zoeNewFeatureVec (void) {
	zoeFeatureVec vec = zoeMalloc(sizeof(struct zoeFeatureVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushFeatureVec (zoeFeatureVec vec, const zoeFeature f) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(zoeFeature));
	}
	vec->elem[vec->size] = zoeCopyFeature(f);
	vec->last = vec->elem[vec->size];
	vec->size++;
}

zoeFeatureVec zoeCopyFeatureVec (const zoeFeatureVec vec) {
	int i;
	zoeFeatureVec v = zoeNewFeatureVec();
	for (i = 0; i < vec->size; i++) zoePushFeatureVec(v, vec->elem[i]);
	return v;
}

#endif
