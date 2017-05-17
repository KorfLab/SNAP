/******************************************************************************\
zoeState.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_STATE_H
#define ZOE_STATE_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zoeFeature.h"
#include "zoeTools.h"

typedef enum {
	INTERNAL,
	EXTERNAL,
	SHUTTLE
} zoeStateType;

struct zoeState  {
	zoeStateType type;      /* enumerated above */
	zoeLabel     label;     /* typical state label */
	float        init;      /* initial probability of state */
	float        term;      /* terminal probability of state */
	int          min;       /* minimum length */
	int          max;       /* maximum length, -1 is unlimited */
	int          geometric; /* true/false for geometric length */
};
typedef struct zoeState * zoeState;

void     zoeDeleteState (zoeState);
zoeState zoeNewState (zoeStateType, zoeLabel, float, float, int, int, int);
zoeState zoeReadState (FILE *);
void     zoeWriteState (FILE *, const zoeState);

#endif
