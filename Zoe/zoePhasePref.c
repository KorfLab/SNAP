/******************************************************************************\
zoePhasePref.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_PHASEPREF_C
#define ZOE_PHASEPREF_C

#include "zoePhasePref.h"

void zoeDeletePhasePref (zoePhasePref pp) {	
	if (pp == NULL) return;
	zoeFree(pp);
	pp = NULL;
}

zoePhasePref zoeNewPhasePref (void) {
	int          i;
	zoePhasePref pp = zoeMalloc(sizeof(struct zoePhasePref));
	
	for (i = 0; i < zoePHASECOUNT; i++) {
		pp->score[i] = MIN_SCORE;
		pp->prob[i]  = 0;
	}
	return pp;
}

zoePhasePref zoeReadPhasePref (FILE * stream) {
	int          i;
	zoePhasePref pp = zoeNewPhasePref();
	
	for (i = 0; i < zoePHASECOUNT; i++) {
		if (fscanf(stream, "%f", &pp->prob[i]) != 1) {
			zoeWarn("zoeReadPhasePref fscanf error");
			zoeDeletePhasePref(pp);
			return NULL;
		}
		/* convert to score too */
		pp->score[i] = zoeFloat2Score(pp->prob[i]);
	}
	return pp;
}

void zoeWritePhasePref (FILE *stream, const zoePhasePref pp) {
	int i;
	
	for (i = 0; i < zoePHASECOUNT; i++) {
		zoeS(stream, "%f\n", pp->prob[i]);
	}
}

score_t zoeScorePhase (zoePhasePref pp, zoeLabel from, zoeLabel to, int inc5) {
	
	switch (from) {
		case Einit:
			switch (to) {
				case Int0:                           return pp->score[Ei_I0];
				case Int1: case Int1T:               return pp->score[Ei_I1];
				case Int2: case Int2TA: case Int2TG: return pp->score[Ei_I2];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		case Exon:
			switch (inc5) {
				case 0: /* corresponds to E0 */
					switch (to) {
						case Int0:                           return pp->score[E0_I0];
						case Int1: case Int1T:               return pp->score[E0_I1];
						case Int2: case Int2TA: case Int2TG: return pp->score[E0_I2];	
						default: zoeExit("zoeScorePhase impossible A %d", to);
					}
				case 1: /* corresponds to E2 */
					switch (to) {
						case Int0:                           return pp->score[E2_I0];
						case Int1: case Int1T:               return pp->score[E2_I1];
						case Int2: case Int2TA: case Int2TG: return pp->score[E2_I2];
						default: zoeExit("zoeScorePhase impossible B %d", to);
					}
				case 2: /* corresponds to E1 */
					switch (to) {
						case Int0:                           return pp->score[E1_I0];
						case Int1: case Int1T:               return pp->score[E1_I1];
						case Int2: case Int2TA: case Int2TG: return pp->score[E1_I2];
						default: zoeExit("zoeScorePhase impossible C %d", to);
					}
				default: zoeExit("zoeScorePhase impossible D");
			}
		case Int0:
			switch (to) {
				case Exon:  return pp->score[I0_E0];
				case Eterm: return pp->score[I0_Et];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		case Int1: case Int1T:
			switch (to) {
				case Exon:  return pp->score[I1_E1];
				case Eterm: return pp->score[I1_Et];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		case Int2: case Int2TA: case Int2TG:
			switch (to) {
				case Exon:  return pp->score[I2_E2];
				case Eterm: return pp->score[I2_Et];
				default: zoeExit("zoeScorePhase does not allow to %d", to);
			}
		default: return 0;
	}
}

#endif
