/******************************************************************************\
 zoeModel.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_MODEL_C
#define ZOE_MODEL_C

#include "zoeModel.h"

/******************************************************************************\
 PRIVATE FUNCTIONS
\******************************************************************************/

zoeModel zoe_ReadModel (FILE * stream, int header_only, score_t pseudocount) {
	char     type[16];
	char     name[16];
	char     data[16];
	float    score;
	int      length, focus, symbols, submodels, mem, i;
	zoeModel model = NULL;
	
	/* parse the model header */
	if (fscanf(stream, "%s %s %d %d %d %d %f",
			name, type, &length, &focus, &symbols, &submodels, &score) != 7) {
		zoeWarn("zoeReadModel header failed");
		return NULL;
	}
	
	/* set model attributes */
	model = zoeNewModel();
	model->length    = length;
	model->focus     = focus;
	model->symbols   = symbols;
	model->submodels = submodels;
	model->score     = score;
	model->name = zoeMalloc(strlen(name) + 1);
	(void)strcpy(model->name, name);
		
	     if (strcmp(type, "WMM") == 0) model->type = WMM;
	else if (strcmp(type, "LUT") == 0) model->type = LUT;
	else if (strcmp(type, "SAM") == 0) model->type = SAM;
	else if (strcmp(type, "SDT") == 0) model->type = SDT;
	else if (strcmp(type, "CDS") == 0) model->type = CDS;
	else if (strcmp(type, "MIX") == 0) model->type = MIX;
	else if (strcmp(type, "TRM") == 0) model->type = TRM;
	else {
		zoeWarn("zoeReadModel sequence model unknown (%s)", type);
		zoeDeleteModel(model);
		return NULL;
	}
	
	/* parse the model body */
	if (model->type == WMM) {
		mem = length * symbols;
		model->data = zoeMalloc(mem * sizeof(score_t));
		if (header_only) {
			for (i = 0; i < mem; i++) {
				model->data[i] = pseudocount;
			}
		} else {
			for (i = 0; i < mem; i++) {
				if (fscanf(stream, "%s", data) != 1) {
					zoeWarn("zoeReadModel wmm fscanf error");
					zoeDeleteModel(model);
					return NULL;
				}
				model->data[i] = zoeText2Score(data);
			}
		}
	} else if (model->type == LUT) {
		mem = zoePOWER[model->symbols][model->length];
		model->data = zoeMalloc(mem * sizeof(score_t));
		if (header_only) {
			for (i = 0; i < mem; i++) {
				model->data[i] = pseudocount;
			}
		} else {
			for (i = 0; i < mem; i++) {
				if (fscanf(stream, "%s", data) != 1) {
					zoeWarn("zoeReadModel lut fscanf error");
					zoeDeleteModel(model);
					return NULL;
				}
				model->data[i] = zoeText2Score(data);
			}
		}
	} else if (model->type == TRM) {
		model->submodel = NULL;
		model->data = NULL;
	} else {
		model->submodel = zoeMalloc(sizeof(struct zoeModel) * model->submodels);
		model->data = NULL;
		for (i = 0; i < model->submodels; i++) {
			if ((model->submodel[i] = zoe_ReadModel(stream, header_only, pseudocount)) == NULL) {
				zoeWarn("zoeReadModel submodel");
				zoeDeleteModel(model);
				return 0;
			}
		}
	}
	
	model->max_left  = zoeModelLengthLeft(model);
	model->max_right = zoeModelLengthRight(model);
	
	return model;
}

static void zoe_WriteModel(FILE * stream, const zoeModel model, int indent, int header) {
	char typestring[16];
	char data[16];
	int  i, j, p;
		
	switch (model->type) {
		case WMM: strcpy(typestring, "WMM"); break;
		case LUT: strcpy(typestring, "LUT"); break;
		case SAM: strcpy(typestring, "SAM"); break;
		case SDT: strcpy(typestring, "SDT"); break;
		case CDS: strcpy(typestring, "CDS"); break;
		case MIX: strcpy(typestring, "MIX"); break;
		case TRM: strcpy(typestring, "TRM"); break;
		default:  strcpy(typestring, "???"); break;
	}
	
	/* print header */
	for (i = 0; i < indent; i++) (void)fputc('\t', stream);
	zoeS(stream, "%s %s %d %d %d %d",
		model->name, typestring, model->length, model->focus, model->symbols,
		model->submodels);
	zoeScore2Text(model->score, data);
	zoeS(stream, " %s\n", data);
	
	/* print body */
	if (model->type == WMM) {
		if (header) return;
		for (i = 0; i < model->length; i++) {
			for (j = 0; j <= indent; j++) (void)fputc('\t', stream); /* padding */
			for (j = 0; j < model->symbols; j++) {
				zoeScore2Text(model->data[ (i * model->symbols) + j], data);
				zoeS(stream, "%s\t", data);
			}
			zoeS(stream, "\n");
		}
	} else if (model->type == LUT) {
		if (header) return;
		p = zoePOWER[model->symbols][model->length];
		for (i = 0; i < p; i++) {
			if (i % model->symbols == 0) { /* padding */
				for (j = 0; j <= indent; j++) (void)fputc('\t', stream);
			}
			zoeScore2Text(model->data[i], data);
			zoeS(stream, "%s\t", data);
			if ( (i+1) % model->symbols == 0) zoeS(stream, "\n");
		}
	} else if (model->type == TRM) {
		/* no body to model */
	} else {
		for (i = 0; i < model->submodels; i++) {
			(void)zoe_WriteModel(stream, model->submodel[i], indent +1, header);
		}
	}
}

static zoeModel zoe_StringToModel (const char * text, float pseudo) {
	FILE * file;
	zoeModel model;
	
	file = tmpfile();
	if (file == NULL) zoeExit("tmpfile failed");
	zoeS(file, "%s\n", text);
	rewind(file);
	model = zoeReadModelHeader(file, pseudo);
	fclose(file);
	return model;
}

/******************************************************************************\
 PUBLIC FUNCTIONS
\******************************************************************************/

void zoeDeleteModel (zoeModel model) {
	int i;
	
	if (model == NULL) return;
	
	for (i = 0; i < model->submodels; i++) {
		zoeDeleteModel(model->submodel[i]);
		model->submodel[i] = NULL;
	}
	if (model->submodel) {
		zoeFree(model->submodel);
		model->submodel = NULL;
	}
	if (model->name) {
		zoeFree(model->name);
		model->name = NULL;
	}
	if (model->data) {
		zoeFree(model->data);
		model->data = NULL;
	}
	zoeFree(model);
	model = NULL;
}

zoeModel zoeNewModel (void) {
	zoeModel model = zoeMalloc(sizeof(struct zoeModel));
	
	model->type      = UNDEFINED_MODEL;
	model->name      = NULL;
	model->length    = UNDEFINED_COOR;
	model->focus     = UNDEFINED_COOR;
	model->symbols   = 0;
	model->submodels = 0;
	model->submodel  = NULL;
	model->score     = MIN_SCORE;
	model->data      = NULL;
	return model;
}

zoeModel zoeReadModel (FILE * stream) {
	return zoe_ReadModel(stream, 0, 0);
}

zoeModel zoeReadModelHeader (FILE * stream, score_t pseudocount) {
	return zoe_ReadModel(stream, 1, pseudocount);
}

void zoeWriteModel (FILE * stream, const zoeModel model) {
	zoe_WriteModel(stream, model, 0, 0); /* indent 0, header + body */
}

void zoeWriteModelHeader (FILE * stream, const zoeModel model) {
	zoe_WriteModel(stream, model, 0, 1); /* indent 0, header true */
}

void zoeAmbiguateModel (zoeModel model, score_t score) {
	int       i, j, index;
	size_t    mem;
	score_t * d5;
	char      string[50];
	
	if (model->symbols == 0 && model->type == TRM) return; /* OK */
	if (model->symbols != 4) zoeExit("zoeAmbiguate error (%d!=4)", model->symbols); 
	model->symbols = 5;
	
	if (model->type == WMM) {
		mem = model->length * 5 * sizeof(score_t);
		d5 = zoeMalloc(mem);
		
		for (i = 0; i < model->length * 5; i++) d5[i] = score;
		
		for (i = 0; i < model->length; i++) {
			/* record the 4-symbol scores in the 5-symbol array */
			for (j = 0; j < 4; j++) {
				d5[(i * 5) + j] = model->data[(i * 4) + j];
			}
		}
		zoeFree(model->data);
		model->data = d5;
		
	} else if (model->type == LUT) {
		mem = zoePOWER[5][model->length] * sizeof(score_t);
		d5 = zoeMalloc(mem);
		
		/* set all values to score to begin */
		for (i = 0; i < zoePOWER[5][model->length]; i++) d5[i] = score;
		
		/* convert all 4-symbol values into the 5-symbol array */
		for (i = 0; i < zoePOWER[4][model->length]; i++) {
			zoeDecToBase(i, 4, string);      /* convert to base4 string */
			index = zoeBaseToDec(5, string); /* convert to base5 string */
			d5[index] = model->data[i];
		}
		zoeFree(model->data);
		model->data = d5;
		
	} else {
		for (i = 0; i < model->submodels; i++) {
			zoeAmbiguateModel(model->submodel[i], score);
		}
	}
}

void zoeDeambiguateModel (zoeModel model) {
	int       i, j, index;
	size_t    mem;
	score_t * d4;
	char      string[50];
	int       containsN;
	
	if (model->symbols == 0 && model->type == TRM) return; /* OK */
	if (model->symbols != 5) {
		zoeExit("zoeDeambiguateModel requires 5-symbol models");
	}
	model->symbols = 4;
		
	if (model->type == WMM) {
		mem = model->length * 4 * sizeof(score_t);
		d4 = zoeMalloc(mem);
		
		for (i = 0; i < model->length; i++) {
			for (j = 0; j < 4; j++) {
				d4[(i * 4) + j] = model->data[(i * 5) + j];
			}
		}
		zoeFree(model->data);
		model->data = d4;
		
	} else if (model->type == LUT) {
		mem = zoePOWER[4][model->length] * sizeof(score_t);
		d4 = zoeMalloc(mem);
				
		/* convert all 5-symbol values into the 4-symbol array */
		for (i = 0; i < zoePOWER[5][model->length]; i++) {
			zoeDecToBase(i, 5, string); /* base5 string */
			containsN = 0;
			for (j = 0; j < strlen(string); j++) {
				if (string[j] == '4') {
					containsN = 1;
					break;
				}
			}
			if (containsN) continue; /* skip N's */
			index = zoeBaseToDec(4, string);
			d4[index] = model->data[i];
		}
		zoeFree(model->data);
		model->data = d4;
	} else {
		for (i = 0; i < model->submodels; i++) {
			zoeDeambiguateModel(model->submodel[i]);
		}
	}
}

zoeModel zoeGetModel (const char * file) {
	FILE     * stream = NULL;
	zoeModel   model  = NULL;
	
	if ((stream = fopen(file, "r")) == NULL)
		zoeExit("zoeGetModel failed to open %s", file);
	if ((model = zoeReadModel(stream)) == NULL)
		zoeExit("zoeGetModel failed to parse %s", file);
	
	(void)fclose(stream);
	return model;
}

zoeModel zoeGetModelHeader (const char * file, score_t pseudocounts) {
	FILE     * stream = NULL;
	zoeModel   model  = NULL;
	
	if ((stream = fopen(file, "r")) == NULL)
		zoeExit("zoeGetModel failed to open %s", file);
	if ((model = zoeReadModelHeader(stream, pseudocounts)) == NULL)
		zoeExit("zoeGetModel failed to parse %s", file);
	
	(void)fclose(stream);
	return model;
}

coor_t zoeModelLengthLeft (const zoeModel model) {
	coor_t length, max = 0;
	int    i;
	
	
	length = model->focus;
	if (length > max) max = length;
	
	if (model->submodels) {	
		for (i = 0; i < model->submodels; i++) {
			length = zoeModelLengthLeft(model->submodel[i]);
			if (length > max) max = length;
		}
	}
	
	return max;
}

coor_t zoeModelLengthRight (const zoeModel model) {
	coor_t length, max = 0;
	int    i;
	
	length = model->length - model->focus -1;
	if (length > max) max = length;
	
	if (model->submodels) {
		for (i = 0; i < model->submodels; i++) {
			length = zoeModelLengthRight(model->submodel[i]);
			if (length > max) max = length;
		}
	}
	
	return max;
}


zoeModel zoeNewCodingModel (int order, float pseudo) {
	char text[256];
	int bytes = 0;

	if (order < 2) zoeExit("coding order must be >= 2");
	bytes += sprintf(text, "Coding CDS 3 2 4 3 0\n");
	bytes += sprintf(text+bytes, "\tframe0 LUT %d %d 4 0 0\n", order +1, order);
	bytes += sprintf(text+bytes, "\tframe1 LUT %d %d 4 0 0\n", order +1, order);
	bytes += sprintf(text+bytes, "\tframe2 LUT %d %d 4 0 0\n", order +1, order);
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewIntronModel (int order, float pseudo) {
	char text[64];
	
	sprintf(text, "Intron LUT %d %d 4 0 0\n", order +1, order);
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewInterModel (int order, float pseudo) {
	char text[64];
	
	sprintf(text, "Inter LUT %d %d 4 0 0\n", order +1, order);
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewAcceptorModel (int order, int length, float pseudo) {
	char text[8192];
	int bytes = 0;
	int i;
	
	bytes += sprintf(text, "Acceptor SDT 2 1 4 2 0\n");
	if (order == 0) {
		bytes += sprintf(text+bytes, "\tAG WMM %d %d 4 0 0\n", length, length -4);
	} else {
		bytes += sprintf(text+bytes, "\tAG SAM %d %d 4 %d 0\n", length, length -4, length);
		for (i = 0; i < length; i++) {
			if      (length -i == 5) bytes += sprintf(text+bytes, "\t\tA LUT %d %d 4 0 0\n", order +1, order);
			else if (length -i == 4) bytes += sprintf(text+bytes, "\t\tG LUT %d %d 4 0 0\n", order +1, order);
			else if (length -i >  3) bytes += sprintf(text+bytes, "\t\tI-%d LUT %d %d 4 0 0\n", length -i -3, order +1, order);
			else                     bytes += sprintf(text+bytes, "\t\tE+%d LUT %d %d 4 0 0\n", i -length +4, order +1, order);
		}
	}
	bytes += sprintf(text+bytes, "\tNN TRM 0 0 0 0 0\n");

	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewDonorModel (int order, int length, float pseudo) {
	char text[8192];
	int bytes = 0;
	int i;
	
	bytes += sprintf(text, "Donor SDT 2 0 4 2 0\n");
	if (order == 0) {
		bytes += sprintf(text+bytes, "\tGT WMM %d 3 4 0 0\n", length);
	} else {
		bytes += sprintf(text+bytes, "\tGT SAM %d 3 4 %d 0\n", length, length);
		for (i = 0; i < length; i++) {
			if      (i == 3) bytes += sprintf(text+bytes, "\t\tG LUT %d %d 4 0 0\n", order +1, order);
			else if (i == 4) bytes += sprintf(text+bytes, "\t\tT LUT %d %d 4 0 0\n", order +1, order);
			else if (i <  3) bytes += sprintf(text+bytes, "\t\tE-%d LUT %d %d 4 0 0\n", 3 -i, order +1, order);
			else             bytes += sprintf(text+bytes, "\t\tI+%d LUT %d %d 4 0 0\n", i -4, order +1, order);
		}
	}
	bytes += sprintf(text+bytes, "\tNN TRM 0 0 0 0 0\n");

		
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewStartModel (int order, int length, float pseudo) {
	char text[8192];
	int bytes = 0;
	int i;
	
	bytes += sprintf(text+bytes, "Start SDT 3 0 4 2 0\n");
	if (order == 0) {
		bytes += sprintf(text+bytes, "\tATG WMM %d %d 4 0 0\n", length, length -6);
	} else {
		bytes += sprintf(text+bytes, "\tATG SAM %d %d 4 %d 0\n", length, length -6, length);
		for (i = 0; i < length; i++) {
			if      (length -i == 6) bytes += sprintf(text+bytes, "\t\tA LUT %d %d 4 0 0\n", order +1, order);
			else if (length -i == 5) bytes += sprintf(text+bytes, "\t\tT LUT %d %d 4 0 0\n", order +1, order);
			else if (length -i == 4) bytes += sprintf(text+bytes, "\t\tG LUT %d %d 4 0 0\n", order +1, order);
			else if (length -i >  3) bytes += sprintf(text+bytes, "\t\tN-%d LUT %d %d 4 0 0\n", length -i -6, order +1, order);
			else                     bytes += sprintf(text+bytes, "\t\tE+%d LUT %d %d 4 0 0\n", i -length +4, order +1, order);
		}
	}
	bytes += sprintf(text+bytes, "\tNNN TRM 0 0 0 0 0\n");
	
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewStopModel (int length, float pseudo) {
	char text[8192];
	int bytes = 0;
	
	bytes += sprintf(text, "Stop SDT 3 0 4 4 0\n");
	bytes += sprintf(text+bytes, "\tTAA WMM %d 6 4 0 0\n", length);
	bytes += sprintf(text+bytes, "\tTAG WMM %d 6 4 0 0\n", length);
	bytes += sprintf(text+bytes, "\tTGA WMM %d 6 4 0 0\n", length);
	bytes += sprintf(text+bytes, "\tNNN TRM 0 0 0 0 0\n");
	
	/*
	bytes += sprintf(text, "Stop SDT 3 0 4 2 0\n");
	bytes += sprintf(text+bytes, "\tTRR SDT 3 0 4 2 0\n");
	bytes += sprintf(text+bytes, "\t\tTGG WMM 1 0 4 0 0\n");
	bytes += sprintf(text+bytes, "\t\tNNN WMM %d 3 4 0 0\n", length);
	bytes += sprintf(text+bytes, "\tNNN TRM 1 0 4 0 0\n");
	*/
	
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewUTR5Model (int order, float pseudo) {
	char text[8192];
	
	sprintf(text, "UTR5 LUT %d %d 4 0 0\n", order +1, order);
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewUTR3Model (int order, float pseudo) {
	char text[8192];
	
	sprintf(text, "UTR3 LUT %d %d 4 0 0\n", order +1, order);
	return zoe_StringToModel(text, pseudo);
}

zoeModel zoeNewPolyAModel (int order, int length, float pseudo) {
	char text[8192];
	int bytes = 0;
	int i;
		
	if (order == 0) {
		bytes += sprintf(text, "PolyA WMM %d 0 4 0 0\n", length);
	} else {
		bytes += sprintf(text, "PolyA SAM %d 0 4 %d 0\n", length, length);
		for (i = 0; i < length; i++) {
			bytes += sprintf(text+bytes, "\t\t%d LUT %d %d 4 0 0\n", i, order +1, order);
		}
	}
	return zoe_StringToModel(text, pseudo); 
}

#endif
