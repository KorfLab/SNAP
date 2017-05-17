/******************************************************************************\
zoePhasePref.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_PHASEPREF_H
#define ZOE_PHASEPREF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zoeMath.h"
#include "zoeFeature.h"
#include "zoeTools.h"

#define zoePHASECOUNT 18

typedef enum {
	Ei_I0, /* transition from Einit to Int0 */
	Ei_I1, /* transition from Einit to Int1 */
	Ei_I2, /* transition from Einit to Int2 */
	E0_I0, /* transition from Exon0 to Int0 */
	E0_I1, /* transition from Exon0 to Int1 */
	E0_I2, /* transition from Exon0 to Int2 */
	E1_I0, /* transition from Exon1 to Int0 */
	E1_I1, /* transition from Exon1 to Int1 */
	E1_I2, /* transition from Exon1 to Int2 */
	E2_I0, /* transition from Exon2 to Int0 */
	E2_I1, /* transition from Exon2 to Int1 */
	E2_I2, /* transition from Exon2 to Int2 */
	
	I0_E0, /* transition from Int0 to Exon0 */
	I0_Et, /* transition from Int0 to Eterm */
	I1_E1, /* transition from Int1 to Exon1 */
	I1_Et, /* transition from Int1 to Eterm */
	I2_E2, /* transition from Int2 to Exon2 */
	I2_Et  /* transition from Int2 to Eterm */
} zoePhasePrefName;

struct zoePhasePref  {
	score_t score[zoePHASECOUNT];
	float   prob[zoePHASECOUNT];
};
typedef struct zoePhasePref * zoePhasePref;

void         zoeDeletePhasePref (zoePhasePref);
zoePhasePref zoeNewPhasePref (void);
zoePhasePref zoeReadPhasePref (FILE *);
void         zoeWritePhasePref (FILE *, const zoePhasePref);
score_t      zoeScorePhase (zoePhasePref, zoeLabel, zoeLabel, int);

#endif
