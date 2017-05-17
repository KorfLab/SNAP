/******************************************************************************\
zoeIsochore.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_ISOCHORE_C
#define ZOE_ISOCHORE_C

#include "zoeHMM.h"
#include "zoeIsochore.h"

void zoeDeleteIsochore (zoeIsochore iso) {
	int i;
	zoeHMM hmm;
	
	if (iso->min_GC) {
		zoeDeleteFVec(iso->min_GC);
		iso->min_GC = NULL;
	}
	if (iso->max_GC) {
		zoeDeleteFVec(iso->max_GC);
		iso->max_GC = NULL;
	}
	if (iso->hmm_file) {
		zoeDeleteTVec(iso->hmm_file);
		iso->hmm_file = NULL;
	}
	if (iso->hmms) {
		for (i = 0; i < iso->hmms->size; i++) {
			hmm = iso->hmms->elem[i];
			zoeDeleteHMM(hmm);
			hmm = NULL;
		}
		zoeDeleteVec(iso->hmms);
	}
}

zoeIsochore zoeNewIsochore (void) {
	zoeIsochore iso = zoeMalloc(sizeof(struct zoeIsochore));
	
	iso->count    = 0;
	iso->min_GC   = zoeNewFVec();
	iso->max_GC   = zoeNewFVec();
	iso->hmm_file = zoeNewTVec();
	iso->hmms     = zoeNewVec();

	return iso;
}


zoeIsochore zoeReadIsochore (FILE * stream) {
	char        name[256];
	int         i, count;
	float       min, max;
	zoeIsochore iso = zoeNewIsochore();
	zoeHMM      hmm;
		
	if (fscanf(stream, "%s %d", name, &count) != 2)
		zoeExit("zoeReadIsochore failed to read header");
	if (strcmp(name, "zoeIsochore") != 0)
		zoeExit("zoeReadIsochore found an unrecognized file type");
		
	iso->count = count;
	for (i = 0; i < count; i++) {
		if (fscanf(stream, "%f %f %s", &min, &max, name) != 3) zoeExit("zoeReadIsochore format error");
		zoePushFVec(iso->min_GC, min);
		zoePushFVec(iso->max_GC, max);
		zoePushTVec(iso->hmm_file, name);
		hmm = zoeGetHMM(name);
		zoePushVec(iso->hmms, hmm);
	}
	
	return iso;
}

zoeIsochore zoeGetIsochore (const char * file) {
	FILE        * stream = NULL;
	zoeIsochore   iso    = NULL;
	char        * ZOE    = getenv("ZOE");
	char          path[1024];
	
	stream = fopen(file, "r");
	if (stream == NULL) {
		sprintf(path, "%s/HMM/%s", ZOE, file);
		stream = fopen(path, "r");
		if (stream == NULL) {
			zoeExit("error opening isochore file");
		}
	}
	
	iso = zoeReadIsochore(stream);
	if (iso == NULL) zoeExit("error reading isochore file");
	
	(void)fclose(stream);
	return(iso);
}

zoeHMM zoeSelectIsochore (const zoeIsochore iso, float GC_fraction) {
	int    i;
	float  min, max;
	zoeHMM hmm;
		
	for (i = 0; i < iso->hmms->size; i++) {
		min = iso->min_GC->elem[i];
		max = iso->max_GC->elem[i];
		hmm = iso->hmms->elem[i];
		if (GC_fraction > min && GC_fraction < max) return hmm;
	}
	
	zoeExit("zoeSelectIsochore could not find a valid isochore for %f", GC_fraction);
	return NULL;
}

#endif
