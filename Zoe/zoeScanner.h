/******************************************************************************\
zoeScanner.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_SCANNER_H
#define ZOE_SCANNER_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zoeDNA.h"
#include "zoeFeature.h"
#include "zoeMath.h"
#include "zoeModel.h"
#include "zoeTools.h"

struct zoeScanner  {
	coor_t               min_pos;    /* minimum scoring position */
	coor_t               max_pos;    /* maximum scoring position */
	zoeDNA               dna;        /* scanners require a seq */
	zoeDNA               anti;       /* reverse-complement */
	zoeModel             model;      /* scanners require a model */
	struct zoeScanner ** subscanner; /* for higher order models */
	char               * sig;        /* binary signature (SDT) */
	score_t            * uscore;     /* user-defined score */
	score_t            * ascore;     /* user-defined anti-parallel score */
	score_t           (* score) (struct zoeScanner *, coor_t);
	score_t           (* scoref)(struct zoeScanner *, zoeFeature);
};
typedef struct zoeScanner * zoeScanner;

void       zoeDeleteScanner (zoeScanner);
zoeScanner zoeNewScanner (zoeDNA, zoeDNA, zoeModel);
void       zoeSetScannerScore (zoeScanner, coor_t, score_t);

#endif
