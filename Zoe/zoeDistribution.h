/******************************************************************************\
zoeDistribution.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_DISTRIBUTION_H
#define ZOE_DISTRIBUTION_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zoeMath.h"
#include "zoeTools.h"

typedef enum {
	UNKNOWN,
	DEFINED,
	GEOMETRIC,
	POISSON,
	CONSTANT
} zoeDistributionType;

struct zoeDistribution  {
	zoeDistributionType   type;   /* the type of function (see above) */
	coor_t                start;  /* starting coordinate */
	coor_t                end;    /* ending coordinate */
	int                   params; /* number of parameters to function */
	float               * param;  /* parameters of function */
};
typedef struct zoeDistribution * zoeDistribution;

void            zoeDeleteDistribution (zoeDistribution);
zoeDistribution zoeNewDistribution (zoeDistributionType, coor_t, coor_t, int, const float *);
zoeDistribution zoeReadDistribution (FILE *);
void            zoeWriteDistribution (FILE *, const zoeDistribution);
score_t         zoeScoreDistribution (const zoeDistribution, coor_t);

#endif
