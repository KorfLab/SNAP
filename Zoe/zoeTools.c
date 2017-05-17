/******************************************************************************\
 zoeTools.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_TOOLS_C
#define ZOE_TOOLS_C

#include "zoeTools.h"

const coor_t   UNDEFINED_COOR = -1;
const frame_t  UNDEFINED_FRAME = -1;
const score_t  MIN_SCORE = -FLT_MAX;
const score_t  MAX_SCORE = FLT_MAX;
const strand_t UNDEFINED_STRAND = 0;


/******************************************************************************\
 Library Information
\******************************************************************************/

static char zoeVersionNumber[] = "2017-03-01";

void zoeLibInfo (void) {
	zoeS(stderr, "ZOE library version %s\n", zoeVersionNumber);
}


/******************************************************************************\
 Program Name
\******************************************************************************/

static char PROGRAM_NAME[256] = "unnamed program";

void zoeSetProgramName (const char * string) {
	(void)strcpy(PROGRAM_NAME, string);
}

char * zoeGetProgramName (void) {
	return PROGRAM_NAME;
}

/******************************************************************************\
 Commandline Processing
\******************************************************************************/

static zoeTVec zoeCL_ARGS = NULL; /* command line arguments */
static zoeHash zoeCL_OPTS = NULL; /* command line options */
static zoeHash zoePR_OPTS = NULL; /* program-defined options */

void zoeSetOption (const char * attr, int flag) {
	flag++; /* 1 = flag, 2 = option with paramter */
	
	if (zoePR_OPTS == NULL) {
		zoePR_OPTS = zoeNewHash();
		zoeCL_OPTS = zoeNewHash();
		zoeCL_ARGS = zoeNewTVec();
	}
	
	if (flag < 1 || flag > 2) zoeExit("zoeSetOption requires 0 = flag, 1 = with argument");
	zoeSetHash(zoePR_OPTS, attr, (void *)(size_t)flag);
}

char * zoeOption (const char * attr) {
	return zoeGetHash(zoeCL_OPTS, attr);
}

void zoeParseOptions (int * argc, char ** argv) {
	int    i;
	char * token = NULL;
		
	/* parse command line */
	for (i = 0; i < *argc; i++) {
		token = argv[i];
		if (token[0] == '-' && strlen(token) > 1) {
			/* option found */
			switch ((size_t)zoeGetHash(zoePR_OPTS, token)) {
				case 0:
					zoeExit("zoeParseOptions: unknown option (%s)", token);
					break;
				case 1:
					zoeSetHash(zoeCL_OPTS, token, token);
					break;
				case 2:
					if (*argc == i+1)
						zoeExit("zoeParseOptions: missing argument for %s", token);
					zoeSetHash(zoeCL_OPTS, token, argv[i+1]);
					i++;
					break;
				default:
					zoeExit("impossible instruction");
			}
			
		} else {
			/* does not look like an option, save it for command line */
			zoePushTVec(zoeCL_ARGS, argv[i]);
		}
	}
		
	/* reset command line */
	*argc = zoeCL_ARGS->size;
	for (i = 0; i < zoeCL_ARGS->size; i++) {
		argv[i] = zoeCL_ARGS->elem[i];
	}
}

/******************************************************************************\
 Typedef I/O
\******************************************************************************/

static void reverse_string (char * s) {
	int c, i, j;

	for (i = 0, j = strlen(s) -1; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

static void itoa (int n, char * s) {
	int i, sign;
		
	if ((sign = n) < 0) n = -n;
	i = 0;
	do {
		s[i++] = n % 10 + '0';
	} while ((n /= 10) > 0);
	if (sign < 0) s[i++] = '-';
	s[i] = '\0';
	reverse_string(s);
}

/* coor_t */
void zoeCoor2Text (coor_t val, char * s) {
	if (val == UNDEFINED_COOR) (void)strcpy(s, "...");
	else                       itoa(val, s);
}

coor_t zoeText2Coor (const char * s) {
	int val;
	
	if (strcmp(s, "...") == 0) return UNDEFINED_COOR;
	if ( (sscanf(s, "%d", &val)) != 1) {
		zoeWarn("zoeText2Coor sscanf failure (%s)", s);
		return UNDEFINED_COOR;
	}
	return (coor_t)val;
}

/* frame_t */
void zoeFrame2Text (frame_t val, char * s) {
	if (val == UNDEFINED_FRAME) (void)strcpy(s, ".");
	else                        itoa(val, s);
}

frame_t zoeText2Frame (const char * s) {
	int val;

	if (strcmp(s, ".") == 0) return UNDEFINED_FRAME;
	if ( (sscanf(s, "%d", &val)) != 1) {
		zoeWarn("zoeText2Frame sscanf failure (%s)", s);
		return UNDEFINED_FRAME;
	}
	return (frame_t)val;
}

/* score_t */
void zoeScore2Text (score_t val, char * s) {
	     if (val <= MIN_SCORE) (void)strcpy(s, ".");
	else if (val >= MAX_SCORE) (void)strcpy(s, "*");
	else                       (void)sprintf(s, "%.3f", val);
}

score_t zoeText2Score (const char * s) {
	float val;
	
	if (strcmp(s, ".") == 0) return MIN_SCORE;
	if (strcmp(s, "*") == 0) return MAX_SCORE;
	
	if ( (sscanf(s, "%f", &val)) != 1) {
		zoeWarn("zoeText2Score sscanf failure (%s)", s);
		return MIN_SCORE;
	}
	return (score_t)val;
}

/* strand_t */
void zoeStrand2Text (strand_t val, char * s) {
	if      (val == '+') (void)strcpy(s, "+");
	else if (val == '-') (void)strcpy(s, "-");
	else                 (void)strcpy(s, ".");
}

strand_t zoeText2Strand (const char * s) {
	if (strcmp(s, ".") == 0) return UNDEFINED_STRAND;
	if (strcmp(s, "+") == 0) return '+';
	if (strcmp(s, "-") == 0) return '-';
	
	zoeWarn("zoeText2Strand sscanf failure (%s)", s);
	return UNDEFINED_STRAND;
}


/******************************************************************************\
 Printing and Error Messages
\******************************************************************************/

void zoeS (FILE * stream, const char * fmt, ...) {
	va_list args;
	
	va_start(args, fmt);
	(void)vfprintf(stream, fmt, args);
	va_end(args);
	(void)fflush(stream);	
}

void zoeO (const char * fmt, ...) {
	va_list args;
		
	va_start(args, fmt);
	(void)vfprintf(stdout, fmt, args);
	va_end(args);
	(void)fflush(stdout);	
}

void zoeE (const char * fmt, ...) {
	va_list args;
		
	va_start(args, fmt);
	(void)vfprintf(stderr, fmt, args);
	va_end(args);	
}

void zoeM (FILE *stream, int argc, ...) {
	int       i;
	va_list   ap;
	char    * s;
		
	va_start(ap, argc);
	for (i = 0; i < argc; i++) {
		s = va_arg(ap, char *);
		fprintf(stream, "%s\n", s);
	}
	va_end(ap);
	
	(void)fflush(stream);
}

void zoeWarn (const char * fmt, ...) {
	va_list args;
		
	(void)fprintf(stderr, "ZOE WARNING (from %s): ", zoeGetProgramName());
	va_start(args, fmt);
	(void)vfprintf(stderr, fmt, args);
	va_end(args);
	(void)fprintf(stderr, "\n");
}

void zoeExit (const char * fmt, ...) {
	va_list args;
		
	(void)fflush(stdout);
	(void)fprintf(stderr, "ZOE ERROR (from %s): ", zoeGetProgramName());
	va_start(args, fmt);
	(void)vfprintf(stderr, fmt, args);
	va_end(args);
	(void)fprintf(stderr, "\n");
	zoeLibInfo();
	exit(2);
}


/******************************************************************************\
 Memory Tools
\******************************************************************************/

void * zoeMalloc (size_t size) {
	void * buffer;
	
	if ((buffer = malloc(size)) == NULL) zoeExit("zoeMalloc");
	return buffer;
}

void * zoeCalloc (size_t nobj, size_t size) {
	void * buffer;

	if ((buffer = calloc(nobj, size)) == NULL) zoeExit("zoeCalloc");
	return buffer;  
}

void * zoeRealloc (void * p, size_t size) {
	void * buffer;
	
	if ((buffer = realloc(p, size)) == NULL) zoeExit("zoeRealloc");
	return buffer;  
}

void zoeFree (void * p) {
	free(p);
	p = NULL;
}


/******************************************************************************\
 Comparison Functions
\******************************************************************************/

int zoeTcmp (const void * a, const void * b) {
	return strcmp( *(char **)a, *(char **)b );
}
int zoeIcmp (const void * a, const void * b) {
	return *(int *)a - *(int *)b;
}
int zoeFcmp (const void * a, const void * b) {
	float f = *(float *)a - *(float *)b;
	     if (f > 0) return  1;
	else if (f < 0) return -1;
	else            return  0;
}


/******************************************************************************\
 Integer Vector
\******************************************************************************/

void zoeDeleteIVec (zoeIVec vec) {	
	if (vec == NULL) return;
	if (vec->elem) {
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeIVec zoeNewIVec (void) {
	zoeIVec vec = zoeMalloc(sizeof(struct zoeIVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushIVec (zoeIVec vec, int val) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(int));
	}
	vec->elem[vec->size] = val;
	vec->last = val;
	vec->size++;
}

/******************************************************************************\
 Float Vector
\******************************************************************************/

void zoeDeleteFVec (zoeFVec vec) {
	if (vec == NULL) return;
	if (vec->elem) {
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeFVec zoeNewFVec (void) {
	zoeFVec vec = zoeMalloc(sizeof(struct zoeFVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushFVec (zoeFVec vec, float val) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(float));
	}
	vec->elem[vec->size] = val;
	vec->last = val;
	vec->size++;
}

/******************************************************************************\
 Text Vector
\******************************************************************************/

void zoeDeleteTVec (zoeTVec vec) {
	int i;
	
	if (vec == NULL) return;
	if (vec->elem) {
		for (i = 0; i < vec->size; i++) zoeFree(vec->elem[i]);
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeTVec zoeNewTVec (void) {
	zoeTVec vec = zoeMalloc(sizeof(struct zoeTVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushTVec (zoeTVec vec, const char * text) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(char *));
	}
	vec->elem[vec->size] = zoeMalloc(strlen(text) +1);
	(void)strcpy(vec->elem[vec->size], text);
	vec->last = vec->elem[vec->size];
	vec->size++;
}


/******************************************************************************\
 zoeVec
\******************************************************************************/

void zoeDeleteVec (zoeVec vec) {	
	if (vec == NULL) return;
	if (vec->elem) {
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeVec zoeNewVec (void) {
	zoeVec vec = zoeMalloc(sizeof(struct zoeVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushVec (zoeVec vec, void * thing) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(void *));
	}
	vec->elem[vec->size] = thing;
	vec->last = vec->elem[vec->size];
	vec->size++;
}

/******************************************************************************\
 Generic Hash

The hash function is my own creature. I used the multiplication method as
described  in Cormen et al. but added a 7 periodic component so that for
example, ATG and TGA don't hash to the same index. Because DNA might be hashed,
I thought it would be a good idea to avoid some multipe of 3 and because number
systems are 2 or 10 based, I stayed away from those too. The 7 values chosen
were kind of arbitrary. I have tested this on some rather large text files, and
the hash function separation is quite good even after several levels of
re-hashing. Performance is good too, about 6x faster than the C++ map and 3
times faster than a Perl hash.

\******************************************************************************/

static double zoeHASH_MULTIPLIER[7] = {
	3.1415926536, /* PI */
	2.7182818285, /* e */
	1.6180339887, /* golden mean */
	1.7320508076, /* square root of 3 */
	2.2360679775, /* square root of 5 */
	2.6457513111, /* square root of 7 */
	3.3166247904, /* square root of 11 */
};

static float zoeMAX_HASH_DEPTH = 2.0;         /* The hash will remain between */

static int zoeHashLevelToSlots (int level) {  /* half full and twice-filled */	
	return pow(4, level);                     /* with these values */
}

static int zoeHashFunc (const zoeHash hash, const char * key) {
	int    i;
	double sum;
	
	sum = 0;
	for (i = 0; i < strlen(key); i++) {
		sum += key[i] * zoeHASH_MULTIPLIER[i % 7];
	}
	
	return (int) (hash->slots * (sum - floor(sum)));
}

static void zoeExpandHash (zoeHash hash) {
	int      i, j;
	char   * key = NULL;
	void   * val = NULL;
	int      oldslots = hash->slots;
	zoeVec * oldkey = hash->key;
	zoeVec * oldval = hash->val;
	zoeVec   kvec;
	zoeVec   vvec;
	zoeTVec  keys;
		
	/* create the new hash */
	hash->level = hash->level +1;
	hash->slots = zoeHashLevelToSlots(hash->level);
	hash->key   = zoeMalloc(hash->slots * sizeof(struct zoeVec));
	hash->val   = zoeMalloc(hash->slots * sizeof(struct zoeVec));
	for (i = 0; i < hash->slots; i++) {
		hash->key[i] = zoeNewVec();
		hash->val[i] = zoeNewVec();
	}
	
	/* brand new hash? */
	if (hash->keys->size == 0) return;

	keys = hash->keys;
	hash->keys = zoeNewTVec();
	
	/* transfer old stuff to new hash */
	for (i = 0; i < oldslots; i++) {
		kvec = oldkey[i];
		vvec = oldval[i];
		for (j = 0; j < kvec->size; j++) {
			key = kvec->elem[j];
			val = vvec->elem[j];
			zoeSetHash(hash, key, val);
		}
	}
	
	/* free old stuff */
	for (i = 0; i < oldslots; i++) {
		kvec = oldkey[i];
		vvec = oldval[i];
		zoeDeleteVec(kvec);
		zoeDeleteVec(vvec);
	}
	zoeFree(oldkey);
	zoeFree(oldval);
	zoeDeleteTVec(keys);
	
}

void zoeDeleteHash (zoeHash hash) {
	int i;
	
	if (hash == NULL) return;
	
	for (i = 0; i < hash->slots; i++) {
		if (hash->key[i]) {
			zoeDeleteVec(hash->key[i]);
			hash->key[i] = NULL;
		}
		if (hash->val[i]) {
			zoeDeleteVec(hash->val[i]);
			hash->val[i] = NULL;
		}
	}
	zoeDeleteTVec(hash->keys);
	hash->keys = NULL;
	zoeDeleteVec(hash->vals);
	hash->vals = NULL;
	zoeFree(hash->key);
	hash->key = NULL;
	zoeFree(hash->val);
	hash->val = NULL;
	zoeFree(hash);
	hash = NULL;
}

zoeHash zoeNewHash (void) {
	zoeHash hash = zoeMalloc(sizeof(struct zoeHash));
	hash->level = 0;
	hash->slots = 0;
	hash->keys  = zoeNewTVec();
	hash->vals  = zoeNewVec();
	hash->key   = NULL;
	hash->val   = NULL;
	zoeExpandHash(hash);
	return hash;
}

void * zoeGetHash (const zoeHash hash, const char * key) {
	int    i, index;
	char * string = NULL;

	index = zoeHashFunc(hash, key);
	/* resolve collisions */
	for (i = 0; i < hash->key[index]->size; i++) {
		string = hash->key[index]->elem[i];
		if (strcmp(key, string) == 0) {
			return hash->val[index]->elem[i];
		}
	}
	return NULL; /* return is NULL if not found */
}

void zoeSetHash (zoeHash hash, const char * key, void * val) {
	int    i, index;
	char * string = NULL;
	int    new_key = 1;
	
	index = zoeHashFunc(hash, key);
	
	/* reassign unless new key */
	for (i = 0; i < hash->key[index]->size; i++) {
		string = hash->key[index]->elem[i];
		if (strcmp(key, string) == 0) {
			hash->val[index]->elem[i] = val;
			new_key = 0;
			return;
		}
	}
	
	if (new_key) {
		zoePushTVec(hash->keys, key);
		zoePushVec(hash->key[index], hash->keys->last);
		zoePushVec(hash->vals, val);
		zoePushVec(hash->val[index], hash->vals->last);
	}
	
	/* check if we have to expand the hash */
	if ((float)hash->keys->size / (float)hash->slots >= zoeMAX_HASH_DEPTH) {
		zoeExpandHash(hash);
	}
}

zoeTVec zoeKeysOfHash (const zoeHash hash) {
	int     i;
	zoeTVec vec = zoeNewTVec();
	
	for (i = 0; i < hash->keys->size; i++) zoePushTVec(vec, hash->keys->elem[i]);
	
	return vec;
}

zoeVec zoeValsOfHash (const zoeHash hash) {
	int    i;
	zoeVec vec = zoeNewVec();

	for (i = 0; i < hash->vals->size; i++) zoePushVec(vec, hash->vals->elem[i]);

	return vec;
}

void zoeStatHash (const zoeHash hash) {
	int i, max, min, total, count;

	max = 0;
	min = 1000000;
	total = 0;
	for (i = 0; i < hash->slots; i++) {
		count = hash->val[i]->size;
		total += count;
		if (count > max) max = count;
		if (count < min) min = count;
	}
	zoeS(stdout, "HashStats: level=%d slots=%d keys=%d min=%d max=%d ave=%f\n",
		hash->level, hash->slots, hash->keys->size, min, max,
		(float)total / (float)hash->slots);
}


/******************************************************************************\
 Suffix Tree
\******************************************************************************/

void zoeDeleteXnode (zoeXnode xn) {
	zoeExit("Xnodes are only freed by Xtree");
}

zoeXnode zoeNewXnode (char c) {
	zoeXnode xn = zoeMalloc(sizeof(struct zoeXnode));
	xn->children = zoeNewVec();
	xn->data     = NULL;
	xn->c        = c;
	return xn;
}
zoeXnode zoeSearchXnode (const zoeXnode xn, char c) {
	int      i;
	zoeXnode child;
	for (i = 0; i < xn->children->size; i++) {
		child = xn->children->elem[i];
		if (child->c == c) {
			return child;
		}
	}
	return NULL;
}

void zoeDeleteXtree (zoeXtree xt) {
	int      i;
	zoeXnode node;
		
	/* deallocate memory from all the nodes */
	for (i = 0; i < xt->alloc->size; i++) {
		node = xt->alloc->elem[i];
		zoeDeleteVec(node->children);
		zoeFree(node);
	}
	zoeDeleteVec(xt->alloc);
	
	zoeFree(xt);
}

zoeXtree zoeNewXtree (void) {
	int i;
	
	zoeXtree xt = zoeMalloc(sizeof(struct zoeXtree));
	for (i = 0; i < 256; i++) {
		xt->head[i] = NULL;
	}
	xt->alloc = zoeNewVec();
	
	return xt;
}

void zoeSetXtree (zoeXtree xt, const char * string, void * value) {
	int      i, len;
	char     c;
	zoeXnode parent, child;
		
	/* string length check */
	len = strlen(string);
	if (len < 1) {
		zoeWarn("zoeSetXtree string length must be at least 1");
		return;
	}
	
	/* head node */
	c = string[0];
	if (xt->head[(int)c] == NULL) {
		/* need to make a new head node */
		xt->head[(int)c] = zoeNewXnode(c);
		zoePushVec(xt->alloc, xt->head[(int)c]);
	}
	
	/* other nodes */
	parent = xt->head[(int)c];
	for (i = 1; i <	len; i++) {
		c = string[i];
		child = zoeSearchXnode(parent, c);
		if (child == NULL) {
			/* create a new node here */
			child = zoeNewXnode(c);
			zoePushVec(parent->children, child);
			zoePushVec(xt->alloc, child);
		}
		
		parent = child;
	}
	
	/* last node - gets data */
	parent->data = value;	
}

void* zoeGetXtree (const zoeXtree xt, const char * string) {
	int      i, len;
	char     c;
	zoeXnode parent, child;
		
	/* string length check */
	len = strlen(string);
	if (len < 1) {
		zoeWarn("zoeGetXtree string length must be at least 1");
		return NULL;
	}
	
	/* head node */
	c = string[0];
	if (xt->head[(int)c] == NULL) return NULL;
	
	/* other nodes */
	parent = xt->head[(int)c];
	for (i = 1; i <	len; i++) {
		c = string[i];
		child = zoeSearchXnode(parent, c);
		if (child == NULL) return NULL;
		parent = child;
	}
	
	/* last node - return data */
	return parent->data;
}

void zoeXtreeInfo (const zoeXtree xt) {
	int i;
	int head_count = 0;
	
	for (i = 0; i < 256; i++) if (xt->head[i] != NULL) head_count++;
	zoeO("head=%d alloc=%d\n", head_count, xt->alloc->size);
}

/******************************************************************************\
 zoeOpenFile zoeCloseFile
\******************************************************************************/

/*

	0	error opening file
	1	standard file
	2	compressed file

*/

static int zoe_file_type (const char * filename) {
	int length = strlen(filename);
	FILE * stream;
		
	stream = fopen(filename, "r");
	if (stream == NULL) return 0; /* failed to open the file */
	
	/* terminates with .gz */
	if (filename[length -3] == '.' &&
		filename[length -2] == 'g' &&
		filename[length -1] == 'z') return 2;
	
	/* terminates with .z */
	if (filename[length -2] == '.' &&
		filename[length -1] == 'z') return 2;
	
	/* terminates with .Z */
	if (filename[length -2] == '.' &&
		filename[length -1] == 'Z') return 2;
	
	return 1;
}

void zoeCloseFile (zoeFile file) {
	switch (file.type) {
		case 0: zoeExit("file already closed"); break;
		case 1: fclose(file.stream); break;
		case 2: pclose(file.stream); break;
		default: zoeExit("odd file type in zoeCloseFile (%d)", file.type);
	}
	
	file.type    = 0;
	file.stream  = NULL;
	file.name[0] = '\0';
}

zoeFile zoeOpenFile (const char * name) {
	zoeFile file;
	char    command[1024];
	
	file.type    = 0;
	file.stream  = NULL;
	file.name[0] = '\0';
	
	/* type */
	file.type = zoe_file_type(name);
	if (file.type == 0) zoeExit("zoeOpenFile failed to open file (%s)", name);
	
	
	/* stream */
	if (file.type == 1) {
		file.stream = fopen(name, "r");
	} else {
		sprintf(command, "gunzip -c %s", name);
		file.stream = popen(command, "r");
	}
	
	/* name */
	strcpy(file.name, name);

	return file;
}

#endif
