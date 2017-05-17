/******************************************************************************\
zoeProtein.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_PROTEIN_H
#define ZOE_PROTEIN_H

#include <stdlib.h>
#include <string.h>

#include "zoeFastaFile.h"
#include "zoeTools.h"

struct zoeProtein  {
	coor_t   length;
	char   * def;
	char   * seq;
	char   * s22; /* 20 amino acids plus X and stop */
};
typedef struct zoeProtein * zoeProtein;

int        zoe_char2aa (const int);
char       zoe_aa2char (const int);
void       zoeDeleteProtein (zoeProtein);
zoeProtein zoeNewProtein (const char *, const char *);
zoeProtein zoeNewTrustedProtein (const char *, const char *);
void       zoeWriteProtein (FILE *, const zoeProtein);
zoeProtein zoeGetProtein (const char *);

#endif
