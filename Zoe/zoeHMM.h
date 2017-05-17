/******************************************************************************\
zoeHMM.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_HMM_H
#define ZOE_HMM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeDuration.h"
#include "zoeState.h"
#include "zoeModel.h"
#include "zoeFeature.h"
#include "zoePhasePref.h"
#include "zoeTransition.h"
#include "zoeTools.h"

struct zoeHMM  {
	/* general attributes */
	char * name;        /* completely arbitrary */
	int    states;      /* number of states in the HMM */
	int    transitions; /* number of transtions in the HMM */
	int    durations;   /* number of duration models */
	int    models;      /* number of sequence models */
	
	/* object storage */
	zoeState      * state;      /* array of zoeState */
	zoeTransition * transition; /* array of zoeTransition */
	zoeDuration   * duration;   /* array of zoeDuration */
	zoeModel      * model;      /* array of zoeModel */
	zoePhasePref    phasepref;  /* exon-intron phase preferences */
	
	/* object mapping */
	zoeDuration dmap[zoeLABELS];            /* durations */
	zoeModel    mmap[zoeLABELS];            /* models */
	zoeState    smap[zoeLABELS];            /* states */
	zoeIVec     jmap[zoeLABELS][zoeLABELS]; /* jump list (reverse arrows) */
	score_t     imap[zoeLABELS];            /* initial probability */
	score_t     kmap[zoeLABELS];            /* terminal probability */
	score_t     xmap[zoeLABELS];            /* geometric extension score */
	score_t     tmap[zoeLABELS][zoeLABELS]; /* transition score */
	coor_t      cmap[zoeLABELS];            /* coordinate adjustments */
};
typedef struct zoeHMM * zoeHMM;

void   zoeSetNscore (zoeLabel, score_t);
void   zoeSetAscore (zoeLabel, score_t);
void   zoeDeleteHMM (zoeHMM);
zoeHMM zoeNewHMM (void);
zoeHMM zoeReadHMM (FILE *);
void   zoeWriteHMM (FILE *, const zoeHMM);
zoeHMM zoeGetHMM (const char *);
score_t zoeGetAscore (int);

#endif
