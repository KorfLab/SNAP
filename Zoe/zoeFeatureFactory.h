/******************************************************************************\
zoeFeatureFactory.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FEATUREFACTORY_H
#define ZOE_FEATUREFACTORY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeDNA.h"
#include "zoeScanner.h"
#include "zoeFeature.h"
#include "zoeFeatureTable.h"
#include "zoeTools.h"

struct zoeFeatureFactory  {
	/* used by all factories */
	zoeFeatureVec (* create)(struct zoeFeatureFactory *, coor_t);
	zoeLabel         type;
	zoeDNA           dna;
	
	/* RFactory and OFactory/XFactory */
	int        length;  /* min length */
	zoeScanner scanner; /* for those factories that use a scanner */
	
	/* OFactory/XFactory only */
	zoeFeatureVec orfs;
	zoeHash       hash;
	score_t       score;
	strand_t      strand;
	
	/* EFactory */
	int       offset;
	score_t * cds[3]; /* cds score in each FRAME */
	score_t * start;  /* score at each valid position */
	score_t * stop;   /* MIN_SCORE at invalid positions */
	score_t * acc;
	score_t * don;
	int * fstop;  /* position of previous stop in this frame */
	/* changed to int 2004-04-06 */
	
};
typedef struct zoeFeatureFactory * zoeFeatureFactory;

void              zoeDeleteFeatureFactory (zoeFeatureFactory);
zoeFeatureFactory zoeNewEFactory (zoeScanner, zoeScanner, zoeScanner, zoeScanner, zoeScanner);
zoeFeatureFactory zoeNewOFactory (zoeScanner, coor_t, score_t, strand_t);
zoeFeatureFactory zoeNewXFactory (zoeScanner, zoeScanner, zoeScanner, zoeScanner, zoeScanner, coor_t, score_t);
zoeFeatureFactory zoeNewRFactory (zoeScanner, coor_t);
zoeFeatureFactory zoeNewSFactory (zoeScanner, zoeLabel);

#endif
