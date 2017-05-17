/******************************************************************************\
zoeTrellis.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_TRELLIS_H
#define ZOE_TRELLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeCDS.h"
#include "zoeDNA.h"
#include "zoeFeatureFactory.h"
#include "zoeHMM.h"
#include "zoeModel.h"
#include "zoeScanner.h"
#include "zoeFeature.h"
#include "zoeFeatureTable.h"
#include "zoeTools.h"

struct zoeTrellis {
	zoeDNA              dna;
	zoeDNA              anti;
	zoeHMM              hmm;
	zoeFeatureVec       xdef;
	score_t             max_score;           /* set at the end */
	score_t             exp_score;           /* expected score of null model */
	zoeFeatureVec       keep;                /* max features */
	zoeIVec             jump;                /* trace-back */
	int                 min_len[zoeLABELS];  /* minimum length (internal & external) */
	int                 max_len[zoeLABELS];  /* maximum explicit length (internal only) */
	zoeScanner          scanner[zoeLABELS];  /* map hmm models to scanners here */
	zoeFeatureFactory   factory[zoeLABELS];  /* external feature factory */
	int                 internal[zoeLABELS]; /* internal states used */	
	int               * trace[zoeLABELS];    /* viterbi trace-back */
	score_t           * score[zoeLABELS];    /* viterbi score */
	zoeFeatureVec       features[zoeLABELS]; /* current features[state_label] */
	score_t          (* ext)(struct zoeTrellis *, coor_t, zoeLabel, zoeFeature);
};
typedef struct zoeTrellis * zoeTrellis;

void       zoeDeleteTrellis (zoeTrellis);
zoeTrellis zoeNewTrellis (zoeDNA, zoeHMM, zoeFeatureVec);
zoeVec     zoePredictGenes (zoeTrellis);
void       zoeScoreCDS(zoeTrellis, zoeCDS, int, int);
void       zoeSetTrellisMeter (int);
void       zoeSetTrellisPadding (int);
char*      zoeGetPartialProtein (zoeTrellis, zoeLabel, zoeFeature);
score_t    zoeScoreExon   (zoeTrellis, zoeFeature, int, int);
score_t    zoeScoreIntron (zoeTrellis, zoeFeature, int);

#endif
