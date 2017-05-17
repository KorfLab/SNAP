/******************************************************************************\
zoeFeature.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FEATURE_H
#define ZOE_FEATURE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zoeTools.h"

#define zoeLABELS 33

typedef enum {
	None,       /* 0 A non-feature? */
	
	Inter,      /* 1 Intergenic */
	Int0,	    /* 2 phase 0 intron */
	Int1,	    /* 3 phase 1 intron */
	Int1T,      /* 4 phase 1 intron, previous exon has overhanging T */
	Int2,	    /* 5 phase 2 intron */
	Int2TA,     /* 6 phase 2 intron, previous exon has overhanging TA */
	Int2TG,     /* 7 phase 2 intron, previous exon has overhanging TG */
	Intron,     /* 8 intron */

	UTR5,       /* 9 5' UTR */
	UTR3,       /* 10 3' UTR */
	
	Esngl,      /* 11 single exon gene */
	Einit,      /* 12 initial exon */
	Eterm,      /* 13 terminal exon */
	Exon,       /* 14 generic or internal exon */
	Coding,     /* 15 generically coding */
	Gene,       /* 16 geneircally gene */
	
	Acceptor,   /* 17 splice acceptor */
	Donor,      /* 18 splice donor */
	Start,      /* 19 ATG */
	Stop,       /* 20 stop codon */
	
	Repeat,     /* 21 Repetitive element */
	CNS,        /* 22 Conserved sequence */
	ORF,        /* 23 Open Reading Frame */
	
	PolyA,      /* 24 poly-A signal */
	Prom,       /* 25 promoter */
	BPS,        /* 26 branch point signal */
	TSS,        /* 27 Trans-splice site */
	
	Misc,       /* 28 miscellaneous feature */
	
	HSP_NN,     /* 29 high-scoring pair: nt-nt (BLASTN) */
	HSP_NA,     /* 20 high-scoring pair: nt-aa (BLASTX) */
	HSP_AN,     /* 31 high-scoring pair: aa-nn (TBLASTN) */
	HSP_AA      /* 32 high-scoring pair: aa-aa (BLASTP, TBLASTX) */
	
} zoeLabel;

struct zoeFeature  {
	zoeLabel    label;
	coor_t      start;
	coor_t      end;
	strand_t    strand;
	score_t     score;
	frame_t     inc5;
	frame_t     inc3;
	frame_t     frame;
	char      * group;
	/*struct zoeFeatureVec * subfeatures;*/
};
typedef struct zoeFeature * zoeFeature;

struct zoeFeatureVec  {
	zoeFeature * elem;   /* array of features */
	int          size;   /* number of features */
	int          limit;  /* number of elements currently allocated */
	zoeFeature   last;
};
typedef struct zoeFeatureVec * zoeFeatureVec;

void       zoeDeleteFeature (zoeFeature);
zoeFeature zoeNewFeature (zoeLabel, coor_t, coor_t, strand_t, score_t, frame_t, frame_t, frame_t, const char */*, zoeFeatureVec*/);
zoeFeature zoeNewTriteFeature (zoeLabel, coor_t, coor_t, const char *);
void       zoeWriteFeature (FILE *, const zoeFeature);
void       zoeWriteDebugFeature (FILE *, const zoeFeature);
void       zoeWriteTriteFeature (FILE *, const zoeFeature);
void       zoeWriteGFF (FILE *, const zoeFeature, const char *, const char *);
zoeFeature zoeReadFeature (FILE *);
zoeFeature zoeReadGFF (FILE*);
zoeFeature zoeCopyFeature (const zoeFeature);
void       zoeAntiFeature (zoeFeature, int);
int        zoeVerifyFeature (const zoeFeature);
int        zoeFeatureCmp (const zoeFeature, const zoeFeature);
int        zoeFeatureCmpPtr (const void *, const void *);
int        zoeFeaturesOverlap (const zoeFeature, const zoeFeature);

zoeLabel zoeText2Label (const char *);
void     zoeLabel2Text (zoeLabel, char *);
void     zoeWriteLabel (FILE *, zoeLabel);

void          zoeDeleteFeatureVec (zoeFeatureVec);
zoeFeatureVec zoeNewFeatureVec (void);
void          zoePushFeatureVec (zoeFeatureVec, const zoeFeature);
zoeFeatureVec zoeCopyFeatureVec (const zoeFeatureVec);


#endif
