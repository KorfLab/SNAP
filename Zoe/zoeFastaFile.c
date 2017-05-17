/******************************************************************************\
 zoeFastaFile.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_FASTA_FILE_C
#define ZOE_FASTA_FILE_C

#include "zoeFastaFile.h"

void zoeDeleteFastaFile(zoeFastaFile entry) {	
	if (entry == NULL) return;
	
	if (entry->def) {
		zoeFree(entry->def);
		entry->def = NULL;
	}
	if (entry->seq) {
		zoeFree(entry->seq);
		entry->seq = NULL;
	}
	zoeFree(entry);
	entry = NULL;
}

zoeFastaFile zoeNewFastaFile (const char * def, const char * seq) {
	zoeFastaFile ff = zoeMalloc(sizeof(struct zoeFastaFile));
	
	ff->def    = zoeMalloc(strlen(def) + 1);
	ff->seq    = zoeMalloc(strlen(seq) + 1);
	ff->length = strlen(seq);
	
	strcpy(ff->def, def);
	strcpy(ff->seq, seq);
	return ff;
}

zoeFastaFile zoeReadFastaFile (FILE * stream) {
    char           c;
	size_t         size;         /* current allocation */
	int            i;
    char         * seq = NULL;   /* sequence */
    char         * def = NULL;   /* definition */
	zoeFastaFile   entry = NULL;
	
	/* initial check for fasta format */
	c = fgetc(stream);
    if (feof(stream)) {
		return NULL;
	}
	if (c != '>') {
		zoeWarn("zoeReadFastaFile > not found");
		return NULL;
	}
	ungetc(c, stream);
		
	/* read the def line */
	size = 256; /* most definitions are small */
	i = 0;
	def = zoeMalloc(size * sizeof(char));
    
    while ((c = fgetc(stream))) {
        if (feof(stream)){
            break;
        }

		if (c == '\n') break;
		def[i] = c;
		i++;
		if (i == size) {
			size *= 2;
			def = zoeRealloc(def, size);
		}
	}
	def[i] = '\0';	
	
	/* read the sequence */
	size = 65536; /* most sequences are large */
	i = 0;
	seq = zoeMalloc(size * sizeof(char));
    
    while ((c = fgetc(stream))) {
        if (feof(stream)){
            break;
        }

		if (c == '>') {
			ungetc(c, stream);
			break; /* next record found */
		}
		if (isspace((int)c)) continue; /* skip spaces */
		seq[i] = c;
		i++;
		
		if (i == size) {
			size *= 2;
			seq = zoeRealloc(seq, size);
		}
	}
	seq[i] = '\0';
    
    entry = zoeNewFastaFile(def+1, seq); /* '>' stripped off definition */
    zoeFree(def);
    zoeFree(seq);
    return entry;
}

static unsigned int zoeFastaLineLength = 50;

void zoeSetFastaLineLength (unsigned int length) {
	zoeFastaLineLength = length;
}

void zoeWriteFastaFile (FILE * stream, const zoeFastaFile entry) {
	coor_t i;

	if (entry->def[0] != '>') zoeS(stream, ">");

	zoeS(stream, "%s", entry->def);
	if (entry->def[strlen(entry->def) -1] != '\n') zoeS(stream, "\n");
	
	for (i=1; i <= entry->length; i++) {
		fputc((entry->seq)[i-1], stream);
		if ((i % zoeFastaLineLength) == 0) zoeS(stream, "\n");
	}
	if ((i % zoeFastaLineLength) != 1) zoeS(stream, "\n");
}

zoeFastaFile zoeGetFastaFile (const char * filename) {
	FILE * stream;
	zoeFastaFile ff;
	
	if ((stream = fopen(filename, "r")) == NULL) {
		zoeExit("error opening file %s", filename);
	}
	
	ff = zoeReadFastaFile(stream);
	if (ff == NULL) zoeExit("crap, trouble with %s", filename);
	
	return ff;
}

#endif
