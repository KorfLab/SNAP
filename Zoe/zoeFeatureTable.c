/******************************************************************************\
zoeFeatureTable.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FEATURETABLE_C
#define ZOE_FEATURETABLE_C

#include <stdio.h>
#include <string.h>

#include "zoeFeatureTable.h"

/******************************************************************************\
 PRIVATE FUNCTIONS
\******************************************************************************/

static int REGION_SIZE = 1000; /* genomic regions broken into 1,000 bp */

static void zoe_IndexLastFeature (zoeFeatureTable ft) {
	char         string[12]; /* ought to be big enough for whole genomes */
	int          index, low, high;
	void       * vec;
	zoeFeature   f = ft->vec->last;
	
	low  = f->start / REGION_SIZE;
	high = f->end / REGION_SIZE;
		
	for (index = low; index <= high; index++) {
		sprintf(string, "%d", index);
		
		vec = zoeGetHash(ft->region, string);
		if (vec == NULL) {
			zoePushTVec(ft->regions, string);
			zoeSetHash(ft->region, ft->regions->last, zoeNewVec());
			vec = zoeGetHash(ft->region, string);
		}
		zoePushVec(vec, f);
	}	
}


/******************************************************************************\
 PUBLIC FUNCTIONS
\******************************************************************************/

void zoeDeleteFeatureTable (zoeFeatureTable ft) {
	int i;
	zoeVec vec;
	zoeTVec keys;
	char * key;
	
	if (ft == NULL) return;
	
	if (ft->def) {
		zoeFree(ft->def);
		ft->def = NULL;
	}
	
	if (ft->vec) {
		zoeDeleteFeatureVec(ft->vec);
		ft->vec = NULL;
	}
	
	if (ft->region) {
		/* reminder: must also delete all the zoeVec objects in HASH */
		keys = zoeKeysOfHash(ft->region);
		for (i = 0; i < keys->size; i++) {
			key = keys->elem[i];
			vec = zoeGetHash(ft->region, key);
			zoeDeleteVec(vec);
		}
		zoeDeleteHash(ft->region);
		zoeDeleteTVec(keys);
		ft->region = NULL;
	}
	
	if (ft->regions) {
		zoeDeleteTVec(ft->regions);
		ft->regions = NULL;
	}
	
	zoeFree(ft);
	ft = NULL;
}

zoeFeatureTable zoeNewFeatureTable (const char * def, const zoeFeatureVec vec) {
	int           i;
	zoeFeatureTable ft = zoeMalloc(sizeof(struct zoeFeatureTable));

	/* copy def and vec */
	ft->def = zoeMalloc(strlen(def) +1);
	(void)strcpy(ft->def, def);

	ft->vec     = zoeNewFeatureVec();
	ft->region  = zoeNewHash();
	ft->regions = zoeNewTVec();
	
	if (vec) for (i = 0; i < vec->size; i++) {
		zoePushFeatureVec(ft->vec, vec->elem[i]);
		zoe_IndexLastFeature(ft);
	}
	
	return ft;
}

zoeFeatureTable zoeReadFeatureTable (FILE * stream) {
	char            c;
	char            * def;
	int             size, i;
	zoeFeature      f  = NULL;
	zoeFeatureTable ft = NULL;
	zoeFeatureVec   fv = NULL;
	
	/* initial check for fasta-ish format */
	c = fgetc(stream);
	if (c == EOF) return NULL;
	if (c != '>') {
		zoeWarn("zoeReadFeatureTable should start with '>'");
		return NULL;
	}
	ungetc(c, stream);
	
	fv = zoeNewFeatureVec();
	
	/* read def line */
	size = 256; /* most definitions are small */
	i = 0;
	def = zoeMalloc(size * sizeof(char));
    
	while ((c = fgetc(stream)) != EOF) {
		if (c == '\n') break;
		def[i] = c;
		i++;
		if (i == size) {
			size *= 2;
			def = zoeRealloc(def, size);
		}
	}
	def[i] = '\0';
	
	
    /* read the features */
	while ((c = fgetc(stream)) != EOF) {
		ungetc(c, stream);
		if (c == '>') {
			break; /* next record found */
		} else {
			if ((f = zoeReadFeature(stream)) == NULL) {
				break; /* none left in file */
			} else {
				zoePushFeatureVec(fv, f);
				zoeDeleteFeature(f);
			}
		}
	}
	
	ft = zoeNewFeatureTable(def+1, fv);
	zoeDeleteFeatureVec(fv);
	zoeFree(def);
	return ft;
}

void zoeWriteFeatureTable (FILE * stream, const zoeFeatureTable ft) {
	int i;

	if (ft->def[0] != '>') zoeS(stream, ">");
	zoeS(stream, "%s", ft->def);
	if (ft->def[strlen(ft->def) -1] != '\n') zoeS(stream, "\n");
	for (i = 0; i < ft->vec->size; i++) {
		zoeWriteFeature(stream, ft->vec->elem[i]);
	}
}

void zoeWriteTriteFeatureTable (FILE * stream, const zoeFeatureTable ft) {
	int i;

	if (ft->def[0] != '>') zoeS(stream, ">");
	zoeS(stream, "%s", ft->def);
	if (ft->def[strlen(ft->def) -1] != '\n') zoeS(stream, "\n");
	for (i = 0; i < ft->vec->size; i++) {
		zoeWriteTriteFeature(stream, ft->vec->elem[i]);
	}
}

void zoeAddFeature (zoeFeatureTable ft, const zoeFeature f) {
	zoePushFeatureVec(ft->vec, f);
	zoe_IndexLastFeature(ft);
}

void zoeAddFeatures (zoeFeatureTable ft, const zoeFeatureVec fv) {
	int i;
		
	for (i = 0; i < fv->size; i++) {
		zoePushFeatureVec(ft->vec, fv->elem[i]);
		zoe_IndexLastFeature(ft);
	}
}

zoeFeatureTable zoeGetFeatureTable (const char * file) {
	FILE            * stream = NULL;
	zoeFeatureTable   ft     = NULL;
	
	if ((stream = fopen(file, "r")) == NULL)
		zoeExit("zoeGetFeatureTable failed to open %s", file);
	if ((ft = zoeReadFeatureTable(stream)) == NULL)
		zoeExit("zoeGetFeatureTable failed to parse %s", file);
	(void)fclose(stream);
	return ft;
}

zoeFeatureTable zoeSelectExons (const zoeFeatureTable ft) {
	int             i;
	zoeFeature      f;
	zoeFeatureTable exons = zoeNewFeatureTable(ft->def, NULL);
	
	for (i = 0; i < ft->vec->size; i++) {
		f = ft->vec->elem[i];
		switch (f->label) {
			case Einit: case Eterm: case Exon: case Esngl:
				zoeAddFeature(exons, f);
				break;
			default:
				break;
		}
	}
	
	return exons;
}

zoeFeatureTable zoeSelectByGroup (const char * def, const zoeFeatureTable ft, const char * group) {
	int             i;
	zoeFeatureVec   vec = zoeNewFeatureVec();
	zoeFeatureTable f2 = NULL;
	
	for (i = 0; i < ft->vec->size; i++) {
		if (group && ft->vec->elem[i]->group) {
			if (strcmp(group, ft->vec->elem[i]->group) == 0) zoePushFeatureVec(vec, ft->vec->elem[i]);
		} else if (group == NULL) {
			zoePushFeatureVec(vec, ft->vec->elem[i]);
		}
	}
	
	f2 = zoeNewFeatureTable(def, vec);
	zoeDeleteFeatureVec(vec);
	return f2;
}

zoeFeatureTable zoeSelectByLabel (const char * def, const zoeFeatureTable ft, zoeLabel label) {
	int             i;
	zoeFeatureVec   vec = zoeNewFeatureVec();
	zoeFeatureTable f2 = NULL;
	
	for (i = 0; i < ft->vec->size; i++) {
		if (label == ft->vec->elem[i]->label) {
			zoePushFeatureVec(vec, ft->vec->elem[i]);
		}
	}
	
	f2 = zoeNewFeatureTable(def, vec);
	zoeDeleteFeatureVec(vec);
	return f2;
}

void zoeAntiFeatureTable(zoeFeatureTable ft, coor_t length) {
	int i;
	
	for (i = 0; i < ft->vec->size; i++) {
		zoeAntiFeature(ft->vec->elem[i], length);
	}
}

zoeTVec zoeFeatureTableGroups (const zoeFeatureTable ft) {
	int     i;
	zoeTVec tvec;
	zoeHash hash = zoeNewHash();
	
	for (i = 0; i < ft->vec->size; i++) {
		if (ft->vec->elem[i]->group) {
			zoeSetHash(hash, ft->vec->elem[i]->group, (void*)1);
		}
	}
	
	tvec = zoeKeysOfHash(hash);
	zoeDeleteHash(hash);
	
	return tvec;
}

zoeVec zoeGetGenes (const zoeFeatureTable ann, const zoeDNA dna) {
	int               i;
	zoeFeatureTable   ft;
	zoeTVec           names;
	char            * name;
	zoeCDS            cds;
	zoeVec            genes = zoeNewVec();
	
	names = zoeFeatureTableGroups(ann);
	for (i = 0; i < names->size; i++) {
		name = names->elem[i];
		ft = zoeSelectByGroup(name, ann, name);
		cds = zoeNewCDS(name, dna, ft->vec);
		zoePushVec(genes, cds);
		zoeDeleteFeatureTable(ft);
	}
		
	if (genes->size > 1) {
		qsort(genes->elem, genes->size, sizeof(zoeCDS), zoeCDScmpptr);
	}
	
	zoeDeleteTVec(names);

	return genes;
}

zoeFeatureVec zoeGetFeaturesNear (const zoeFeatureTable ft, coor_t start, coor_t end) {
	char          string[12];
	int           i, index, low, high;
	zoeVec        vec;
	zoeFeatureVec store = zoeNewFeatureVec();
	zoeFeatureVec keep;
	zoeFeature    f, last_feature;
	
	low  = start / REGION_SIZE;
	high = end / REGION_SIZE;
		
	for (index = low; index <= high; index++) {
		sprintf(string, "%d", index);
		vec = zoeGetHash(ft->region, string);
		if (vec) {
			for (i = 0; i < vec->size; i++) {
				f = vec->elem[i];
				if (f->end < start || f->start > end) continue;
				zoePushFeatureVec(store, f);
			}
		}
	}
	
	if (store->size <= 1) {
		return store;
	}
	
	/* must prune redundancies - sort, remove dups */
	qsort(store->elem, store->size, sizeof(zoeFeature), zoeFeatureCmpPtr);
	keep = zoeNewFeatureVec();
	last_feature = store->elem[0];
	zoePushFeatureVec(keep, store->elem[0]);
	for (i = 1; i < store->size; i++) {
		if (zoeFeatureCmp(store->elem[i], last_feature) != 0) {
			zoePushFeatureVec(keep, store->elem[i]);
			last_feature = store->elem[i];
		}
	}
	
	zoeDeleteFeatureVec(store);
	return keep;
}

void zoePadFeatureTable(zoeFeatureTable ft, int padding) {
 	int i;
 	
 	for (i = 0; i < ft->vec->size; i++) {
 		ft->vec->elem[i]->start += padding;
 		ft->vec->elem[i]->end   += padding;
 	}
 }

#endif
