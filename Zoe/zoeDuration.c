/******************************************************************************\
zoeDuration.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_DURATION_C
#define ZOE_DURATION_C

#include "zoeDuration.h"

void zoeDeleteDuration(zoeDuration dm) {
	int i;
		
	if (dm == NULL) return;
	
	for (i = 0; i < dm->distributions; i++) {
		if (dm->distribution[i]) {
			zoeDeleteDistribution(dm->distribution[i]);
			dm->distribution[i] = NULL;
		}
	}
	if (dm->distribution) {
		zoeFree(dm->distribution);
		dm->distribution = NULL;
	}
	zoeFree(dm);
	dm = NULL;
}

/*

This is an exceptional case! The durations copy only the pointers to the
distributions. But they free the distributions!!

*/

zoeDuration zoeNewDuration (zoeLabel label, int ds, const zoeDistribution * d) {
	int         i;
	zoeDuration dm = zoeMalloc(sizeof(struct zoeDuration));
	
	dm->label         = label;
	dm->distributions = ds;
	dm->distribution  = zoeMalloc(ds * sizeof(zoeDistribution));
	for (i = 0; i < ds; i++) {
		dm->distribution[i] = d[i];
	}
	return dm;
}

zoeDuration zoeReadDuration (FILE * stream) {
	char              duration_name[16];
	int               distributions, j;
	zoeLabel         label;
	zoeDistribution * distribution = NULL;


	if (fscanf(stream, "%s %d", duration_name, &distributions) != 2) {
		zoeWarn("fscanf zoeReadDuration header");
		return NULL;
	}
		
	/* label */
	label = zoeText2Label(duration_name);
	if (label == None) {
		zoeWarn("zoeReadDuration name not allowed (%s)", duration_name);
		return NULL;
	}
	
	/* distributions */
	distribution = zoeMalloc(distributions * sizeof(zoeDistribution));
	for (j = 0; j < distributions; j++) {
		if((distribution[j] = zoeReadDistribution(stream)) == NULL) {
			zoeWarn("zoeReadDuration distribution");
			return NULL;
		}
	}
	
	return zoeNewDuration(label, distributions, distribution);
	
}

void zoeWriteDuration (FILE * stream, const zoeDuration dm) {
	int  j;
	char label[16];

	zoeLabel2Text(dm->label, label);
	zoeS(stream, "%s %d\n", label, dm->distributions);
	for (j = 0; j < dm->distributions; j++) {
		zoeWriteDistribution(stream, dm->distribution[j]);
	}
}

score_t zoeScoreDuration(const zoeDuration dm, coor_t pos) {
	int i, found;

	if (pos < 1) {
		zoeWarn("zoeScoreDuration positive integers only");
		return MIN_SCORE;
	}
	
	found = -1;
	for (i = 0; i < dm->distributions; i++) {
		if (pos >= dm->distribution[i]->start && pos <= dm->distribution[i]->end) {
			found = i;
			break;
		} else if (pos >= dm->distribution[i]->start && dm->distribution[i]->end == 0) {
			found = i;
			break;
		} else if (pos <= dm->distribution[i]->end && dm->distribution[i]->start == 0) {
			found = i;
			break;
		}
	}
	
	if (found == -1) {
		zoeWriteDuration(stderr, dm);
		zoeExit("zoeScoreDuration out of bounds %d", pos);
	}
	return zoeScoreDistribution(dm->distribution[found], pos);
}

#endif
