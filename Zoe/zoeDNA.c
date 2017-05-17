/******************************************************************************\
 zoeDNA.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_DNA_C
#define ZOE_DNA_C

#include "zoeDNA.h"

static char s5CodonTable[5][5][5] = {
	{
		{'K','N','K','N','X',},
		{'T','T','T','T','T',},
		{'R','S','R','S','X',},
		{'I','I','M','I','X',},
		{'X','X','X','X','X',},
	},
	{
		{'Q','H','Q','H','X',},
		{'P','P','P','P','P',},
		{'R','R','R','R','R',},
		{'L','L','L','L','L',},
		{'X','X','X','X','X',},
	},
	{
		{'E','D','E','D','X',},
		{'A','A','A','A','A',},
		{'G','G','G','G','G',},
		{'V','V','V','V','V',},
		{'X','X','X','X','X',},
	},
	{
		{'*','Y','*','Y','X',},
		{'S','S','S','S','S',},
		{'*','C','W','C','X',},
		{'L','F','L','F','X',},
		{'X','X','X','X','X',},
	},
	{
		{'X','X','X','X','X',},
		{'X','X','X','X','X',},
		{'X','X','X','X','X',},
		{'X','X','X','X','X',},
		{'X','X','X','X','X',},
	},

};

static void zoe_s5_stats (zoeDNA dna) {
	int i;
	
	for (i = 0; i < 5; i++) dna->f5[i] = 0;
	for (i = 0; i < 5; i++) dna->c5[i] = 0;
	if (dna->length == 0) return;
	for (i = 0; i < dna->length; i++) dna->c5[(int)dna->s5[i]]++;
	for (i = 0; i < 5; i++) dna->f5[i] = dna->c5[i]/dna->length;
}

char* zoeTranslateS5 (const char *s5, int tx_len, frame_t offset) {
	int          i, idx;
	coor_t       aa_length;
	char       * seq;
	
	if (offset < 0 || offset > 2) {
		zoeExit("attempt to use illegal offset (%d) in zoeTranslateS5\n", offset);
	}

	aa_length = (tx_len - offset) / 3;
	seq = zoeMalloc((size_t)aa_length +1);
	for (i = 0; i < aa_length; i++) {
		idx    = 3 * i + offset;
		seq[i] = s5CodonTable[(int)s5[idx]][(int)s5[idx+1]][(int)s5[idx+2]];
	}
	seq[aa_length] = '\0';
	return seq;
}

void zoeDeleteDNA(zoeDNA dna) {		
	if (dna == NULL) return;
	if (dna->def) {zoeFree(dna->def); dna->def = NULL;}
	if (dna->seq) {zoeFree(dna->seq); dna->seq = NULL;}
	if (dna->s5)  {zoeFree(dna->s5);  dna->s5  = NULL;}
	if (dna->s16) {zoeFree(dna->s16); dna->s16 = NULL;}
	
	zoeFree(dna);
	dna = NULL;
}

zoeDNA zoeNewDNA (const char * def, const char * seq) {
	coor_t i;
	zoeDNA dna = zoeMalloc(sizeof(struct zoeDNA));
	
	dna->length = strlen(seq);	
 	dna->seq    = zoeMalloc(dna->length +1);
	dna->s5     = zoeMalloc(dna->length +1);
	dna->s16    = zoeMalloc(dna->length +1);
	dna->def    = zoeMalloc(strlen(def) +1);
	
	/* set sequence and definition */
	strcpy(dna->def, def);
	strcpy(dna->seq, seq);
	
	/* create s5 sequence and warn on alphabet errors */
	for (i = 0; i < dna->length; i++) {
		switch (seq[i]) {
			case 'A': case 'a': dna->s5[i] = 0; break;
			case 'C': case 'c': dna->s5[i] = 1; break;
			case 'G': case 'g': dna->s5[i] = 2; break;
			case 'T': case 't': dna->s5[i] = 3; break;
			case 'R': case 'r': case 'Y': case 'y':
			case 'M': case 'm': case 'K': case 'k':
			case 'W': case 'w': case 'S': case 's':
			case 'B': case 'b': case 'D': case 'd':
			case 'H': case 'h': case 'V': case 'v':
			case 'N': case 'n': dna->s5[i] = 4; break;
			case '-':           dna->s5[i] = 4; break;
			default:
				dna->s5[i]  = 4;
				dna->seq[i] = 'N';
				zoeWarn("zoeNewDNA editing illegal symbol '%c' to 'N' in %s at %d/%d",
					seq[i], def, i, dna->length);
		}
	}

	/* create s16 sequence but don't warn on alphabet errors */
	for (i = 0; i < dna->length; i++) {
		switch (seq[i]) {
			case 'A': case 'a':	dna->s16[i] =  8; break; /* 1000 */
			case 'C': case 'c':	dna->s16[i] =  4; break; /* 0100 */
			case 'G': case 'g':	dna->s16[i] =  2; break; /* 0010 */
			case 'T': case 't':	dna->s16[i] =  1; break; /* 0001 */
			case 'R': case 'r':	dna->s16[i] = 10; break; /* 1010 */
			case 'Y': case 'y':	dna->s16[i] =  5; break; /* 0101 */
			case 'M': case 'm':	dna->s16[i] = 12; break; /* 1100 */
			case 'K': case 'k':	dna->s16[i] =  3; break; /* 0011 */
			case 'W': case 'w':	dna->s16[i] =  9; break; /* 1001 */
			case 'S': case 's':	dna->s16[i] =  6; break; /* 0110 */
			case 'B': case 'b':	dna->s16[i] =  7; break; /* 0111 */
			case 'D': case 'd':	dna->s16[i] = 11; break; /* 1011 */
			case 'H': case 'h':	dna->s16[i] = 13; break; /* 1101 */
			case 'V': case 'v':	dna->s16[i] = 14; break; /* 1110 */
			case 'N': case 'n':	dna->s16[i] = 15; break; /* 1111 */
			case '-':           dna->s16[i] = 15; break;
			default:            dna->s16[i] = 15;
		}
	}
	
	zoe_s5_stats(dna);

	return dna;
}

zoeDNA zoeCopyDNA(const zoeDNA dna) {
	/* should be re-written to copy memory rather than go through constructor */
	return zoeNewDNA(dna->def, dna->seq);
}

zoeDNA zoeReverseDNA (const char * def, const zoeDNA dna) {
	coor_t   i;
	char   * seq = zoeMalloc(dna->length +1);
	zoeDNA   rev = NULL;

	for (i = 0; i < dna->length; i++) {
		seq[dna->length -i -1] = dna->seq[i];
	}
	seq[dna->length] = '\0';
	rev = zoeNewDNA(def, seq);
	zoeFree(seq);
	return rev;
}

zoeDNA zoeComplementDNA (const char * def, const zoeDNA dna) {
	coor_t   i;
	char   * seq = zoeMalloc(dna->length + 1);
	zoeDNA   comp = NULL;
	
	for (i = 0; i < dna->length; i++) {
		switch (dna->seq[i]) {
			case 'A': seq[i] = 'T'; break;
			case 'C': seq[i] = 'G'; break;
			case 'G': seq[i] = 'C'; break;
			case 'T': seq[i] = 'A'; break;
			case 'R': seq[i] = 'Y'; break;
			case 'Y': seq[i] = 'R'; break;
			case 'M': seq[i] = 'K'; break;
			case 'K': seq[i] = 'M'; break;
			case 'W': seq[i] = 'W'; break;
			case 'S': seq[i] = 'S'; break;
			case 'B': seq[i] = 'V'; break;
			case 'D': seq[i] = 'H'; break;
			case 'H': seq[i] = 'D'; break;
			case 'V': seq[i] = 'B'; break;
			case 'N': seq[i] = 'N'; break;
			case 'a': seq[i] = 't'; break;
			case 'c': seq[i] = 'g'; break;
			case 'g': seq[i] = 'c'; break;
			case 't': seq[i] = 'a'; break;
			case 'r': seq[i] = 'y'; break;
			case 'y': seq[i] = 'r'; break;
			case 'm': seq[i] = 'k'; break;
			case 'k': seq[i] = 'm'; break;
			case 'w': seq[i] = 'w'; break;
			case 's': seq[i] = 's'; break;
			case 'b': seq[i] = 'v'; break;
			case 'd': seq[i] = 'h'; break;
			case 'h': seq[i] = 'd'; break;
			case 'v': seq[i] = 'b'; break;
			case 'n': seq[i] = 'n'; break;
			case '-': seq[i] = '-'; break;
			default:
				zoeExit("hard to reach error in zoeComplementDNA");
		}
	}
	seq[dna->length] = '\0';
	comp = zoeNewDNA(def, seq);
	zoeFree(seq);
	return comp;
}

zoeDNA zoeAntiDNA (const char * def, const zoeDNA dna) {
	zoeDNA rev  = zoeReverseDNA(def, dna);
	zoeDNA anti = zoeComplementDNA(def, rev);
	zoeDeleteDNA(rev);
	return anti;
}

void zoeLCmask (zoeDNA dna) {
	int i;
	
	/* change lowercase to N in s5 and s16 */
	for (i = 0; i < dna->length; i++) {
		if (islower((int)dna->seq[i])) {
			dna->s5[i] = 4;
			if (dna->s16) dna->s16[i] = 15;
		}
	}
	
	zoe_s5_stats(dna);
}

void zoeLCunmask (zoeDNA dna) {
	int i;
	
	for (i = 0; i < dna->length; i++) {
		if (dna->s5[i] == 4) {
			switch (dna->seq[i]) {
				case 'a': dna->s5[i] = 0; break;
				case 'c': dna->s5[i] = 1; break;
				case 'g': dna->s5[i] = 2; break;
				case 't': dna->s5[i] = 3; break;
				default:  dna->s5[i] = 4; break;
			}
		}
	}
	
	zoe_s5_stats(dna);
}

void zoeLCfilter (zoeDNA dna) {
	int i;
	
	/* change lowercase to N */
	for (i = 0; i < dna->length; i++) {
		if (islower((int)dna->seq[i])) {
			dna->s5[i] = 4;
			dna->seq[i] = 'N';
			if (dna->s16) dna->s16[i] = 15;
		}
	}
	
	zoe_s5_stats(dna);
}

void zoeLCsmooth (zoeDNA dna, coor_t flank, coor_t island, coor_t min_len) {
	int           i, j, d, len1, len2;
	zoeFeature    misc, f, prev; 
	zoeFeatureVec fvec;
	char          * s5;
	int           ncount = 0;
	
	/*
		Purpose: Mask long-ish lowercase regions and fill in small
		islands of uppercase. Short repeats may be over-zealous masking.
	*/
	
	
	/* create a copy of the dna->s5 to work with */
	s5 = zoeMalloc(dna->length +1);	
	for (i = 0; i < dna->length; i++) {
		if (dna->s5[i] == 4 || islower((int)dna->seq[i])) {
			s5[i] = 4;
			ncount++;
		} else {
			s5[i] = dna->s5[i];
		}
	}
	
	/* abort if no N's/lowercase in sequence */
	if (ncount == 0) {
		zoeFree(s5);
		return;
	}
	
	/* find all regions of N/lowercase */
	misc = zoeNewFeature(Misc, 0, 0, '+', 0, 0, 0, 0, NULL/*, NULL*/);
	fvec = zoeNewFeatureVec();
	
	for (i = 0; i < dna->length; i++) if (s5[i] == 4) break;
	for (i = 0; i < dna->length; i++) {
		if (s5[i] == 4) {
			for (j = i+1; j < dna->length; j++) {
				if (s5[j] != 4) break;
			}
			misc->start = i;
			misc->end   = j -1;
			zoePushFeatureVec(fvec, misc);
			i = j;
		}
	}
	if (fvec->size == 0) zoeExit("zoeLCmask vec->size 0");
		
	/*
		convert small islands to N
		
		   prev                f
		NNNNNNNNNN acgtgaa NNNNNNNNNN
		<- len1 ->    d    <- len2 ->
	*/
	
	for (i = 1; i < fvec->size; i++) {
		f    = fvec->elem[i];
		prev = fvec->elem[i-1];
		len1 = prev->end - prev->start + 1;
		len2 = f->end - f->start + 1;
		d = f->start - prev->end -1;
			
		if (len1 >= flank && len2 >= flank && d <= island) {
			for (j = prev->end+1; j <= f->start -1; j++) {
				s5[j]  = 4;
			}
		}	
	}
	
	/* clear and start over */
	zoeDeleteFeatureVec(fvec);
	fvec = zoeNewFeatureVec();
	for (i = 0; i < dna->length; i++) if (s5[i] == 4) break;
	for (i = 0; i < dna->length; i++) {
		if (s5[i] == 4) {
			for (j = i+1; j < dna->length; j++) {
				if (s5[j] != 4) break;
			}
			misc->start = i;
			misc->end   = j -1;
			if (j - i > min_len) zoePushFeatureVec(fvec, misc);
			i = j;
		}
	}
	
	/* edit DNA */
	for (i = 0; i < fvec->size; i++) {
		f = fvec->elem[i];
		for (j = f->start; j <= f->end; j++) {
			dna->seq[j] = 'N';
			dna->s5[j] = 4;
			if (dna->s16) dna->s16[j] = 15;
		}
	}
	
	zoeDeleteFeature(misc);
	zoeDeleteFeatureVec(fvec);
	
	zoe_s5_stats(dna);
}

zoeDNA zoeSubseqDNA (const char *def, const zoeDNA dna, coor_t from, coor_t length) {
	zoeDNA   sub = NULL;
	char   * seq = NULL;
	int      i;
		
	seq = zoeMalloc(length +1);
	for (i = from; i < from + length; i++) {
		seq[i - from] = dna->seq[i];
	}
	seq[length] = '\0';	
	sub = zoeNewDNA(def, seq);

	zoeFree(seq);
	return sub;
}

zoeDNA zoeFeatureDNA (const char *def, const zoeDNA dna, const zoeFeature f) {
	zoeDNA sub  = NULL;
	zoeDNA anti = NULL;
	
	sub = zoeSubseqDNA(def, dna, f->start, f->end - f->start +1);
	if (f->strand == '-') {
		anti = zoeAntiDNA(def, sub);
		zoeDeleteDNA(sub);
		return anti;
	} else {
		return sub;
	}
}


zoeProtein zoeTranslateDNA (const char *def, const zoeDNA dna, frame_t offset) {
	char       * seq;
	zoeProtein   pro = NULL;
		
	if (offset < 0 || offset > 2) {
		zoeExit("attempt to use illegal offset (%d) in zoeTranslateDNA in %s\n", offset, def);
	}
	
	seq = zoeTranslateS5(dna->s5, dna->length, offset);
	pro = zoeNewTrustedProtein(def, seq);
	zoeFree(seq);
	return pro;
}

zoeProtein zoeTranslateFeature (const char * def, const zoeDNA dna, const zoeFeature f) {
	zoeDNA     sub;
	zoeProtein pro;
	
	sub = zoeFeatureDNA(def, dna, f);
	pro = zoeTranslateDNA(def, sub, f->inc5);
	zoeDeleteDNA(sub);
	return pro;
}


void zoeWriteDNA (FILE * stream, const zoeDNA dna) {
	coor_t i;

	if (dna->def[0] != '>') zoeS(stream, ">");
	zoeS(stream, "%s", dna->def);
	if (dna->def[strlen(dna->def) -1] != '\n') zoeS(stream, "\n");
	for (i = 1; i <= dna->length; i++) {
		(void)fputc(dna->seq[i-1], stream);
		if ((i % 60) == 0) zoeS(stream, "\n");
	}
	if ((i % 60) != 1) zoeS(stream, "\n");
}

zoeDNA zoeGetDNA (const char * file) {
	FILE         * stream = NULL;
	zoeFastaFile   fasta  = NULL;
	zoeDNA         dna    = NULL;
	
	if ((stream = fopen(file, "r")) == NULL)
		zoeExit("zoeGetDNA failed to open %s", file);
	if ((fasta = zoeReadFastaFile(stream)) == NULL)
		zoeExit("zoeGetDNA failed to parse %s", file);
	(void)fclose(stream);
	
	dna = zoeNewDNA(fasta->def, fasta->seq);
	zoeDeleteFastaFile(fasta);
	
	return(dna);
}

zoeFeatureVec zoeORFs (const zoeDNA dna, strand_t strand) {
	int             i, frame, length;
	zoeDNA          anti   = NULL;
	zoeFeatureVec   orfs   = NULL;
	zoeProtein      pro[3];
	zoeIVec         stops  = NULL;
	zoeFeature      orf;
	coor_t          start, end;
	
	orfs = zoeNewFeatureVec();
	
	/* plus-strand ORFs */
	if (strand == '+' || strand == '=') {
		for (frame = 0; frame < 3; frame++) {
			if ((pro[frame] = zoeTranslateDNA(">zoeORFs translate", dna, frame)) == NULL) {
				zoeExit("zoeORFs fatal error");
			}
			stops = zoeNewIVec();
			for (i = 0; i < pro[frame]->length; i++) {
				if      (pro[frame]->seq[i] == '*') zoePushIVec(stops, i * 3 + frame);
				else if (pro[frame]->seq[i] == 'X') zoePushIVec(stops, i * 3 + frame);
			}
			for (i = 1; i < stops->size; i++) {
				start = stops->elem[i-1] +3; /* one after stop codon */
				end   = stops->elem[i] -1;   /* one before stop codon */
				length = end - start +1;
				if (length > 0) {
					orf = zoeNewFeature(Coding, start, end, '+', MIN_SCORE, 0, 0, frame, NULL/*, NULL*/);
					zoePushFeatureVec(orfs, orf);
					zoeDeleteFeature(orf);
				}
			}
			zoeDeleteProtein(pro[frame]);
			zoeDeleteIVec(stops);
		}
	}
		
	/* minus-strand ORFs */
	if (strand == '-' || strand == '=') {
		anti = zoeAntiDNA(">zoeORFs anti", dna);
		for (frame = 0; frame < 3; frame++) {
			pro[frame] = zoeTranslateDNA(">zoeORFs translate", anti, frame);
			stops = zoeNewIVec();
			for (i = 0; i < pro[frame]->length; i++) {
				if      (pro[frame]->seq[i] == '*') zoePushIVec(stops, i * 3 + frame);
				else if (pro[frame]->seq[i] == 'X') zoePushIVec(stops, i * 3 + frame);
			}
			for (i = 1; i < stops->size; i++) {
				start  = anti->length -stops->elem[i];
				end    = anti->length -stops->elem[i -1] -4;
				length = end - start +1;
				if (length > 0) {
					orf = zoeNewFeature(Coding, start, end, '-', MIN_SCORE, 0, 0, frame, NULL/*, NULL*/);
					zoePushFeatureVec(orfs, orf);
					zoeDeleteFeature(orf);
				}
			}
			zoeDeleteProtein(pro[frame]);
			zoeDeleteIVec(stops);
		}
		zoeDeleteDNA(anti);
	}
	
	return orfs;
}

void zoeWriteFeatureDNA(FILE *stream, const zoeFeature f, const zoeDNA dna, coor_t extra) {
	zoeFeature copy;
	zoeDNA     sub;
	
	copy = zoeCopyFeature(f);
	copy->start -= extra;
	copy->end   += extra;
	sub = zoeFeatureDNA(">", dna, copy);
	zoeS(stream, "%s\n", sub->seq);
	zoeDeleteDNA(sub);
	zoeDeleteFeature(copy);
}

zoeDNA zoeMakePaddedDNA (const zoeDNA real_dna, int padding) {
	coor_t   i, new_length;
	zoeDNA   dna;
	char   * seq;
		
	new_length = real_dna->length + padding * 2;
	seq = zoeMalloc(new_length +1);
	
	for (i = 0; i < padding; i++) seq[i] = 'N';
	for (i = 0; i < real_dna->length; i++) seq[i+padding] = real_dna->seq[i];
	for (i = 0; i < padding; i++) seq[i + real_dna->length + padding] = 'N';
	seq[new_length] = '\0';
	dna = zoeNewDNA(real_dna->def, seq);
	
	zoeFree(seq);
	return dna;
}

#endif
