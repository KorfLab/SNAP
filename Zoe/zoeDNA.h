/******************************************************************************\
zoeDNA.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_DNA_H
#define ZOE_DNA_H

#include <stdio.h>
#include <string.h>

#include "zoeFastaFile.h"
#include "zoeProtein.h"
#include "zoeFeature.h"
#include "zoeTools.h"

/******************************************************************************\

Num  s5    s16
===  ===   ===
 0    A    reserved but not used
 1    C    T  (0001)
 2    G    G  (0010)
 3    T    K  (0011)
 4    N    C  (0100)
 5         Y  (0101)
 6         S  (0110)
 7         B  (0111)
 8         A  (1000)
 9         W  (1001)
10         R  (1010)
11         D  (1011)
12         M  (1100)
13         H  (1101)
14         V  (1110)
15         N  (1111)

\******************************************************************************/

struct zoeDNA  {
	coor_t   length;
	char   * def;     /* definition */
	char   * seq;     /* ascii sequence */
	char   * s5;      /*  5 symbol numeric sequence */
	char   * s16;     /* 15 symbol numeric sequence */
	float    c5[5];   /* symbol counts */
	float    f5[5];   /* symbol frequencies */
};
typedef struct zoeDNA * zoeDNA;

void          zoeDeleteDNA (zoeDNA);
zoeDNA        zoeNewDNA (const char *, const char *);
zoeDNA        zoeCopyDNA (const zoeDNA);
zoeDNA        zoeReverseDNA (const char *, const zoeDNA);
zoeDNA        zoeComplementDNA (const char *, const zoeDNA);
zoeDNA        zoeAntiDNA (const char *, const zoeDNA);
void          zoeLCmask (zoeDNA);
void          zoeLCunmask (zoeDNA);
void          zoeLCfilter (zoeDNA);
void          zoeLCsmooth (zoeDNA, coor_t, coor_t, coor_t);
zoeDNA        zoeSubseqDNA (const char *, const zoeDNA, coor_t, coor_t);
zoeDNA        zoeFeatureDNA (const char *, const zoeDNA, const zoeFeature);
zoeProtein    zoeTranslateDNA (const char *, const zoeDNA, frame_t);
zoeProtein    zoeTranslateFeature (const char *, const zoeDNA, const zoeFeature);
void          zoeWriteDNA (FILE *, const zoeDNA);
zoeDNA        zoeGetDNA (const char *);
zoeFeatureVec zoeORFs (const zoeDNA, strand_t);
void          zoeWriteFeatureDNA(FILE *, const zoeFeature, const zoeDNA, coor_t);
zoeDNA        zoeMakePaddedDNA (const zoeDNA, int);
char*         zoeTranslateS5 (const char*, int, frame_t);

#endif
