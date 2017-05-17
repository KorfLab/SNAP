/******************************************************************************\
 zoeTransition.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_TRANSITION_C
#define ZOE_TRANSITION_C

#include "zoeTransition.h"

void zoeDeleteTransition (zoeTransition t) {
	if (t == NULL) return;
	zoeFree(t);
	t = NULL;
}

zoeTransition zoeNewTransition (const char * from, const char * to, float prob) {
	zoeTransition t = zoeMalloc(sizeof(struct zoeTransition));
	
	t->from  = zoeText2Label(from);
	t->to    = zoeText2Label(to);
	t->prob  = prob;
	t->score = zoeFloat2Score(prob);
	return t;
}

zoeTransition zoeReadTransition (FILE * stream) {
	char  from[16];
	char  to[16];
	float prob;
		
	if (fscanf(stream, "%s %s %f", from, to, &prob) != 3) {
		zoeWarn("fscanf error in zoeReadTransition");
		return NULL;
	}
	return zoeNewTransition(from, to, prob);
}

void zoeWriteTransition (FILE * stream, const zoeTransition t) {
	char from[16];
	char to[16];
		
	zoeLabel2Text(t->from, from);
	zoeLabel2Text(t->to,   to);
	zoeS(stream, "%s\t%s\t%f\n", from, to, t->prob);
}

#endif
