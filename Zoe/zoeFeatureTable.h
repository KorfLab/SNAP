/******************************************************************************\
zoeFeatureTable.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FEATURETABLE_H
#define ZOE_FEATURETABLE_H

#include <stdio.h>
#include <string.h>

#include "zoeCDS.h"
#include "zoeDNA.h"
#include "zoeFeature.h"
#include "zoeTools.h"

struct zoeFeatureTable  {
	char          * def;
	zoeFeatureVec   vec;     /* storage for the actual features */
	zoeHash         region;  /* feature indexing by large-ish region */
	zoeTVec         regions; /* names of indicies */
};
typedef struct zoeFeatureTable * zoeFeatureTable;

void            zoeDeleteFeatureTable (zoeFeatureTable);
zoeFeatureTable zoeNewFeatureTable (const char *, const zoeFeatureVec);
zoeFeatureTable zoeReadFeatureTable (FILE *);
void            zoeWriteFeatureTable (FILE *, const zoeFeatureTable);
void            zoeWriteTriteFeatureTable (FILE *, const zoeFeatureTable);
void            zoeAddFeature (zoeFeatureTable, const zoeFeature);
void            zoeAddFeatures (zoeFeatureTable, const zoeFeatureVec);
zoeFeatureTable zoeSelectExons (const zoeFeatureTable);
zoeFeatureTable zoeSelectByGroup (const char *, const zoeFeatureTable, const char *);
zoeFeatureTable zoeSelectByLabel (const char *, const zoeFeatureTable, zoeLabel);
void            zoeAntiFeatureTable (zoeFeatureTable, coor_t length);
zoeFeatureTable zoeGetFeatureTable (const char *);
zoeTVec         zoeFeatureTableGroups (const zoeFeatureTable);
zoeVec          zoeGetGenes (const zoeFeatureTable, const zoeDNA); /* vector of zoeCDS */
zoeFeatureVec   zoeGetFeaturesNear(const zoeFeatureTable, coor_t, coor_t);
void            zoePadFeatureTable(zoeFeatureTable, int);

#endif

