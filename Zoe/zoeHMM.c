/******************************************************************************\
zoeHMM.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_HMM_C
#define ZOE_HMM_C

#include "zoeHMM.h"

static score_t Nscore[zoeLABELS] = {
	0,                   /* None */
	1,                   /* Inter */
	1, 1, 1, 1, 1, 1, 1, /* Intron and variants */
	-1, -1,                /* UTR5, UTR3 */
	-1, -1, -1, -1, -1,      /* Esngl, Einit, Eterm, Exon, Coding */
	-1, -1, -1, -1,          /* Acceptor, Donor, Start, Stop */
	1, -1, -1,            /* Repeat, CNS, ORF */
	-1, -1, -1, -1,          /* PolyA, Prom, BPS, TSS */
	0,                   /* Misc */
	0, 0, 0, 0           /* HSPs */
};

static score_t Ascore[zoeLABELS] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

score_t zoeGetAscore (int label) {
	return Ascore[label];
}

static int tag_error (FILE * stream, const char * tag) {
	char tmp[32];
	
	if (fscanf(stream, "%s", tmp) != 1) {
		zoeWarn("tag_error %s", tag);
		return 1;
	}
	if (strcmp(tmp, tag) != 0) {
		zoeWarn("tag not found %s", tag);
		return 1;
	}
	return 0;
}

static void zoeMapHMM (zoeHMM hmm) {
	int      i, j, k;
	zoeLabel label, from, to;
	score_t  score;
	score_t  null_score;
	
	/* smap (state map) */
	for (i = 0; i < zoeLABELS; i++) hmm->smap[i] = NULL;
	for (i = 0; i < hmm->states; i++) {
		label = hmm->state[i]->label;
		if (label == Intron) {
			hmm->smap[Int0]   = hmm->state[i];
			hmm->smap[Int1]   = hmm->state[i];
			hmm->smap[Int1T]  = hmm->state[i];
			hmm->smap[Int2]   = hmm->state[i];
			hmm->smap[Int2TA] = hmm->state[i];
			hmm->smap[Int2TG] = hmm->state[i];
		} else {
			hmm->smap[label] = hmm->state[i];
		}
	}

	/* imap (initial score) */
	for (i = 0; i < zoeLABELS; i++) hmm->imap[i] = MIN_SCORE;
	for (i = 0; i < hmm->states; i++) {
		if (hmm->state[i]->type != INTERNAL) continue;
		label = hmm->state[i]->label;
		
		if (hmm->state[i]->init == 0) score = MIN_SCORE;
		else                          score = zoeFloat2Score(hmm->state[i]->init);
		
		if (label == Intron) {
			hmm->imap[Int0]   = score;
			hmm->imap[Int1]   = score;
			hmm->imap[Int1T]  = score;
			hmm->imap[Int2]   = score;
			hmm->imap[Int2TA] = score;
			hmm->imap[Int2TG] = score;
		} else {
			hmm->imap[label] = score;
		}
	}
	
	/* kmap (terminal score) */
	for (i = 0; i < zoeLABELS; i++) hmm->kmap[i] = MIN_SCORE;
	for (i = 0; i < hmm->states; i++) {
		if (hmm->state[i]->type != INTERNAL) continue;
		label = hmm->state[i]->label;
		
		if (hmm->state[i]->term == 0) score = MIN_SCORE;
		else                          score = zoeFloat2Score(hmm->state[i]->term);
		
		if (label == Intron) {
			hmm->kmap[Int0]   = score;
			hmm->kmap[Int1]   = score;
			hmm->kmap[Int1T]  = score;
			hmm->kmap[Int2]   = score;
			hmm->kmap[Int2TA] = score;
			hmm->kmap[Int2TG] = score;
		} else {
			hmm->kmap[label] = score;
		}
	}
	
	/* dmap - duration map */
	for (i = 0; i < zoeLABELS; i++) hmm->dmap[i] = NULL;
	for (i = 0; i < hmm->durations; i++) {
		label = hmm->duration[i]->label;
		if (label == Intron) {
			hmm->dmap[Int0]   = hmm->duration[i];
			hmm->dmap[Int1]   = hmm->duration[i];
			hmm->dmap[Int1T]  = hmm->duration[i];
			hmm->dmap[Int2]   = hmm->duration[i];
			hmm->dmap[Int2TA] = hmm->duration[i];
			hmm->dmap[Int2TG] = hmm->duration[i];
		} else {
			hmm->dmap[label] = hmm->duration[i];
		}
	}
	
	/* xmap (extend score) */
	null_score = zoeScoreDuration(hmm->duration[Inter], 2) - zoeScoreDuration(hmm->duration[Inter], 1);
	for (i = 0; i < zoeLABELS; i++) hmm->xmap[i] = 0;
	for (i = 0; i < hmm->durations; i++) {
		label = hmm->duration[i]->label;
		if (hmm->duration[i]->distribution[0]->type != GEOMETRIC) {
			score = 0; /* explicit duration does not use xmap */
		} else {
			if (hmm->duration[i]->distributions != 1)
				zoeExit("HMM xmap distributions != 1");
			score = zoeScoreDuration(hmm->duration[i], 2) - zoeScoreDuration(hmm->duration[i], 1);
			score -= null_score;
		}
		
		if (label == Intron) {
			hmm->xmap[Int0]   = score;
			hmm->xmap[Int1]   = score;
			hmm->xmap[Int1T]  = score;
			hmm->xmap[Int2]   = score;
			hmm->xmap[Int2TA] = score;
			hmm->xmap[Int2TG] = score;
		} else {
			hmm->xmap[label] = score;
		}
	}
	
	/* mmap - model map */
	for (i = 0; i < zoeLABELS; i++) hmm->mmap[i] = NULL;
	for (i = 0; i < hmm->models; i++) {
		label = zoeText2Label(hmm->model[i]->name);
		if (label == Intron) {
			hmm->mmap[Int0]   = hmm->model[i];
			hmm->mmap[Int1]   = hmm->model[i];
			hmm->mmap[Int1T]  = hmm->model[i];
			hmm->mmap[Int2]   = hmm->model[i];
			hmm->mmap[Int2TA] = hmm->model[i];
			hmm->mmap[Int2TG] = hmm->model[i];
		} else {
			hmm->mmap[label] = hmm->model[i];
		}
	}
	
	/* tmap */
	for (i = 0; i < zoeLABELS; i++) {
		for (j = 0; j < zoeLABELS; j++) hmm->tmap[i][j] = MIN_SCORE;
	}
	for (i = 0; i < hmm->transitions; i++) {
		from  = hmm->transition[i]->from;
		to    = hmm->transition[i]->to;
		score = hmm->transition[i]->score;
		
		if (from == Intron) {
			hmm->tmap[Int0][to]   = score;
			hmm->tmap[Int1][to]   = score;
			hmm->tmap[Int1T][to]  = score;
			hmm->tmap[Int2][to]   = score;
			hmm->tmap[Int2TA][to] = score;
			hmm->tmap[Int2TG][to] = score;
		} else if (to == Intron) {
			hmm->tmap[from][Int0]   = score;
			hmm->tmap[from][Int1]   = score;
			hmm->tmap[from][Int1T]  = score;
			hmm->tmap[from][Int2]   = score;
			hmm->tmap[from][Int2TA] = score;
			hmm->tmap[from][Int2TG] = score;
		} else {
			hmm->tmap[from][to] = score;
		}
	}
	
	/* jmap */
	for (i = 0; i < zoeLABELS; i++) {
		for (j = 0; j < zoeLABELS; j++) {
			hmm->jmap[i][j] = 0;
		}
	}
	for (i = 0; i < zoeLABELS; i++) {
		for (j = 0; j < zoeLABELS; j++) {
			if (hmm->smap[i] == NULL) continue;
			if (hmm->smap[i]->type != INTERNAL) continue;
			for (k = 0; k < zoeLABELS; k++) {
				if (hmm->tmap[i][j] == MIN_SCORE) continue;
				if (hmm->tmap[j][k] == MIN_SCORE) continue;
				if (hmm->jmap[k][j] == NULL) {
					hmm->jmap[k][j] = zoeNewIVec();
				}
				zoePushIVec(hmm->jmap[k][j], i);
			}
		}
	}
	
	/* cmap */
	for (i = 0; i < zoeLABELS; i++) {
		if (hmm->smap[i] == NULL) continue;
		switch (i) {
			case Einit:
				hmm->cmap[i] = hmm->mmap[Start]->max_left + hmm->mmap[Donor]->max_right;
				break;
			case Exon:
				hmm->cmap[i] =  hmm->mmap[Acceptor]->max_left + hmm->mmap[Donor]->max_right;
				break;
			case Eterm:
				hmm->cmap[i] = hmm->mmap[Acceptor]->max_left + hmm->mmap[Stop]->max_right;
				break;
			case Esngl:
				hmm->cmap[i] = hmm->mmap[Start]->max_left + hmm->mmap[Stop]->max_right;
				break;
			case PolyA:
				hmm->cmap[i] = hmm->mmap[PolyA]->max_left + hmm->mmap[PolyA]->max_right;
				break;
			case Prom:
				hmm->cmap[i] = hmm->mmap[Prom]->max_left + hmm->mmap[Prom]->max_right;
				break;
			case TSS:
				hmm->cmap[i] = hmm->mmap[TSS]->max_left + hmm->mmap[TSS]->max_right;
				break;
			case Repeat:
				hmm->cmap[i] = 0;
				break;
			case ORF:
				hmm->cmap[i] = 0;
				break;
			case Int0: case Int1: case Int1T: case Int2: case Int2TA: case Int2TG:
				hmm->cmap[i] = hmm->mmap[Donor]->max_right;
				break;
			case Inter: case UTR5:  case UTR3:
				break;
			default: zoeExit("%d not supported in hmm->cmap", i);
		}
	}

}

static zoeHMM zoeAbortHMM (zoeHMM hmm, const char * string) {
	zoeWarn("zoeHMM aborting (%s)", string);
	zoeDeleteHMM(hmm);
	return NULL;
}

/****************************************************************************\
 PUBLIC FUNCTIONS
\****************************************************************************/

void zoeSetNscore (zoeLabel label, score_t score) {
	Nscore[label] = score;
}

void zoeSetAscore (zoeLabel label, score_t score) {
	Ascore[label] = score;
}

void zoeDeleteHMM (zoeHMM hmm) {
	int i, j;
	
	if (hmm == NULL) return;

	if (hmm->name) {
		zoeFree(hmm->name);
		hmm->name = NULL;
	}
	
	for (i = 0; i < hmm->states; i++) {
		zoeDeleteState(hmm->state[i]);
		hmm->state[i] = NULL;
	}
	if (hmm->state) {
		zoeFree(hmm->state);
		hmm->state = NULL;
	}
	
	for (i = 0; i < hmm->transitions; i++) {
		zoeDeleteTransition(hmm->transition[i]);
		hmm->transition[i] = NULL;
	}
	if (hmm->transition) {
		zoeFree(hmm->transition);
		hmm->transition = NULL;
	}
	
	for (i = 0; i < hmm->durations; i++) {
		zoeDeleteDuration(hmm->duration[i]);
		hmm->duration[i] = NULL;
	}
	if (hmm->duration) {
		zoeFree(hmm->duration);
		hmm->duration = NULL;
	}
	
	for (i = 0; i < hmm->models; i++) {
		zoeDeleteModel(hmm->model[i]);
		hmm->model[i] = NULL;
	}
	if (hmm->model) {
		zoeFree(hmm->model);
		hmm->model = NULL;
	}
			
	/* jmap only mapped needed to free */
	for (i = 0; i < zoeLABELS; i++) {
		for (j = 0; j < zoeLABELS; j++) {
			if (hmm->jmap[i][j]) {
				zoeDeleteIVec(hmm->jmap[i][j]);
			}
		}
	}
	
	zoeFree(hmm);
}

zoeHMM zoeNewHMM (void) {
	int    i, k;
	zoeHMM hmm = zoeMalloc(sizeof(struct zoeHMM));
	
	hmm->name        = NULL;
	hmm->states      = 0;
	hmm->transitions = 0;
	hmm->durations   = 0;
	hmm->models      = 0;
	hmm->state       = NULL;
	hmm->transition  = NULL;
	hmm->duration    = NULL;
	hmm->model       = NULL;
	hmm->phasepref   = NULL;
	
	for (i = 0; i < zoeLABELS; i++) {
		hmm->dmap[i] = NULL;
		hmm->mmap[i] = NULL;
		hmm->smap[i] = NULL;
		hmm->imap[i] = MIN_SCORE;
		hmm->kmap[i] = MIN_SCORE;
		hmm->xmap[i] = MIN_SCORE;
		hmm->cmap[i] = 0;
		
		for (k = 0; k < zoeLABELS; k++) {
			hmm->jmap[i][k] = NULL;
			hmm->tmap[i][k] = MIN_SCORE;
		}		
	}

	return hmm;
}

zoeHMM zoeReadHMM (FILE * stream) {
	char     type[16];
	char     name[256];
	zoeLabel label;
	int      i, states, transitions, durations, models;
	zoeHMM   hmm = zoeNewHMM();
	
	/* header */
	if (fscanf(stream, "%s %s %d %d %d %d", type, name, &states,
			&transitions, &durations, &models) != 6) return zoeAbortHMM(hmm, "header");
	
	/* set attributes */
	if (strcmp(type, "zoeHMM") != 0) return zoeAbortHMM(hmm, "zoeHMM");
	hmm->states      = states;
	hmm->transitions = transitions;
	hmm->durations   = durations;
	hmm->models      = models;
	hmm->name        = zoeMalloc(strlen(name) + 1);
	                   (void)strcpy(hmm->name, name);
	hmm->state      = zoeMalloc(states * sizeof(zoeState));
	hmm->transition = zoeMalloc(transitions * sizeof(zoeTransition));
	hmm->model      = zoeMalloc(models * sizeof(zoeModel));
	hmm->duration   = zoeMalloc(durations * sizeof(zoeDuration));
		
	/* parse the states */
	if (tag_error(stream, "<STATES>"))
		return zoeAbortHMM(hmm, "<STATES>");
	for (i = 0; i < states; i++) {
		if ((hmm->state[i] = zoeReadState(stream)) == NULL)
			return zoeAbortHMM(hmm, "states");
	}
	
	/* parse state transitions */
	if (tag_error(stream, "<STATE_TRANSITIONS>"))
		return zoeAbortHMM(hmm, "<STATE_TRANSITIONS>");
	for (i = 0; i < hmm->transitions; i++) {
		if ((hmm->transition[i] = zoeReadTransition(stream)) == NULL)
			return zoeAbortHMM(hmm, "transitions");
	}
	
	/* parse phase preferences */
	if (tag_error(stream, "<PHASE_PREFERENCES>"))
		return zoeAbortHMM(hmm, "<PHASE_PREFERENCES>");
	if ((hmm->phasepref = zoeReadPhasePref(stream)) == NULL)
		return zoeAbortHMM(hmm, "phasepref");
	
	/* parse state durations */
	if (tag_error(stream, "<STATE_DURATIONS>"))
		return zoeAbortHMM(hmm, "<STATE_DURATIONS>");
	for (i = 0; i < hmm->durations; i++) {
		if ((hmm->duration[i] = zoeReadDuration(stream)) == NULL)
			return zoeAbortHMM(hmm, "durations");
	}
	
	/* parse sequence models */
	if (tag_error(stream, "<SEQUENCE_MODELS>"))
		return zoeAbortHMM(hmm, "<SEQUENCE_MODELS>");
	for (i = 0; i < hmm->models; i++) {
		if ((hmm->model[i] = zoeReadModel(stream)) == NULL) {
			zoeE("model %d failure\n", i);
			return zoeAbortHMM(hmm, "models");
		}
		
		label = zoeText2Label(hmm->model[i]->name);
		if (hmm->model[i]->symbols == 4) {
			zoeAmbiguateModel(hmm->model[i], Nscore[label]);
		}
	}
	zoeMapHMM(hmm);
	return hmm;
}

void zoeWriteHMM (FILE * stream, const zoeHMM hmm) {
	int i;

	zoeS(stream, "zoeHMM\t%s\t%d\t%d\t%d\t%d\n", hmm->name, hmm->states,
		hmm->transitions, hmm->durations, hmm->models);
	zoeS(stream, "\n<STATES>\n");
	for (i = 0; i < hmm->states; i++) {
		zoeWriteState(stream, hmm->state[i]);
	}
	
	zoeS(stream, "\n<STATE_TRANSITIONS>\n");
	for (i = 0; i < hmm->transitions; i++) {
		zoeWriteTransition(stream, hmm->transition[i]);
	}
	
	zoeS(stream, "\n<PHASE_PREFERENCES>\n");
	zoeWritePhasePref(stream, hmm->phasepref);
	
	zoeS(stream, "\n<STATE_DURATIONS>\n");
	for (i = 0; i < hmm->durations; i++) {
		zoeWriteDuration(stream, hmm->duration[i]);
	}

	zoeS(stream, "\n<SEQUENCE_MODELS>\n");
	for (i = 0; i < hmm->models; i++) {
		zoeWriteModel(stream, hmm->model[i]);
	}
}

zoeHMM zoeGetHMM (const char * file) {
	FILE   * stream = NULL;
	zoeHMM   hmm    = NULL;
	char   * ZOE    = getenv("ZOE");
	char     path[1024];
	
	stream = fopen(file, "r");
	if (stream == NULL) {
		sprintf(path, "%s/HMM/%s", ZOE, file);
		stream = fopen(path, "r");
		if (stream == NULL) {
			zoeExit("error opening HMM file (%s)", path);
		}
	}
	hmm = zoeReadHMM(stream);
	if (hmm == NULL) zoeExit("zoeGetHMM failed to parse %s", file);
	
	(void)fclose(stream);
	return(hmm);
}

#endif
