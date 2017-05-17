/******************************************************************************\
 zoeProtein.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_PROTEIN_C
#define ZOE_PROTEIN_C

#include "zoeProtein.h"

int zoe_char2aa (const int c) {
	switch (c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'D': return 2;
		case 'E': return 3;
		case 'F': return 4;
		case 'G': return 5;
		case 'H': return 6;
		case 'I': return 7;
		case 'K': return 8;
		case 'L': return 9;
		case 'M': return 10;
		case 'N': return 11;
		case 'P': return 12;
		case 'Q': return 13;
		case 'R': return 14;
		case 'S': return 15;
		case 'T': return 16;
		case 'V': return 17;
		case 'W': return 18;
		case 'Y': return 19;
		case '*': return 20;
		case 'X': return 21;
		default:  return -1;
	}
}

char zoe_aa2char (const int aa) {
	switch (aa) {
		case  0: return 'A';
		case  1: return 'C';
		case  2: return 'D';
		case  3: return 'E';
		case  4: return 'F';
		case  5: return 'G';
		case  6: return 'H';
		case  7: return 'I';
		case  8: return 'K';
		case  9: return 'L';
		case 10: return 'M';
		case 11: return 'N';
		case 12: return 'P';
		case 13: return 'Q';
		case 14: return 'R';
		case 15: return 'S';
		case 16: return 'T';
		case 17: return 'V';
		case 18: return 'W';
		case 19: return 'Y';
		case 20: return '*';
		case 21: return 'X';
		default: return -1;
	}
}

void zoeDeleteProtein(zoeProtein pro) {	
	if (pro == NULL) return;
	if (pro->def) {
		zoeFree(pro->def);
		pro->def = NULL;
	}
	if (pro->seq) {
		zoeFree(pro->seq);
		pro->seq = NULL;
	}
	if (pro->s22) {
		zoeFree(pro->s22);
		pro->s22 = NULL;
	}
	zoeFree(pro);
	pro = NULL;
}

zoeProtein zoeNewTrustedProtein (const char * def, const char * seq) {
	int  i;
	char aa;
	zoeProtein pro = zoeMalloc(sizeof(struct zoeProtein));
	
	pro->length = strlen(seq);
	pro->def    = zoeMalloc(strlen(def) +1);
	pro->seq    = zoeMalloc(pro->length +1);
	pro->s22    = zoeMalloc(pro->length +1);
	strcpy(pro->def, def);
	strcpy(pro->seq, seq);
	for (i = 0; i < strlen(seq); i++) {
		aa = zoe_char2aa(seq[i]);
		pro->s22[i] = aa;
	}
	return pro;
}

zoeProtein zoeNewProtein (const char * def, const char * seq) {
	int        i;
	char       aa;
	zoeProtein pro = zoeMalloc(sizeof(struct zoeProtein));
	
	pro->length = strlen(seq);
	pro->def    = zoeMalloc(strlen(def) +1);
	pro->seq    = zoeMalloc(pro->length +1);
	pro->s22    = zoeMalloc(pro->length +1);
	strcpy(pro->def, def);
	strcpy(pro->seq, seq);
	
	for (i = 0; i < strlen(seq); i++) {
		aa = zoe_char2aa(seq[i]);
		if (aa == -1) zoeWarn("illegal amino acid (%c)", seq[i]);
		pro->s22[i] = aa;
	}
	
	
	return pro;
}

void zoeWriteProtein (FILE * stream, const zoeProtein pro) {
	coor_t i;
	
	if (pro->def[0] != '>') zoeS(stream, ">");
	zoeS(stream, "%s", pro->def);
	if (pro->def[strlen(pro->def) -1] != '\n') zoeS(stream, "\n");

	for (i = 1; i <= pro->length; i++) {
		(void)fputc(pro->seq[i-1], stream);
		if ((i % 50) == 0) zoeS(stream, "\n");
	}
	if ((i % 50) != 1) zoeS(stream, "\n");
}

zoeProtein zoeGetProtein (const char * file) {
	FILE         * stream = NULL;
	zoeFastaFile   fasta  = NULL;
	zoeProtein     pro    = NULL;
	
	if ((stream = fopen(file, "r")) == NULL)
		zoeExit("zoeGetProtein failed to open %s", file);
	if ((fasta = zoeReadFastaFile(stream)) == NULL)
		zoeExit("zoeGetProtein failed to parse %s", file);
	pro = zoeNewProtein(fasta->def, fasta->seq);
	
	zoeDeleteFastaFile(fasta);
	(void)fclose(stream);
	
	return(pro);
}


#endif
