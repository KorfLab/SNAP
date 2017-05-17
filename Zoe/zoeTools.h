/******************************************************************************\
 zoeTools.h - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_TOOLS_H
#define ZOE_TOOLS_H

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char * zoeFunction;
const char * zoeConstructor;
const char * zoeMethod;

void   zoeLibInfo (void);
void   zoeSetProgramName (const char*);
char * zoeGetProgramName (void);
void   zoeSetOption(const char *, int);
void   zoeParseOptions (int *, char **);
char * zoeOption (const char *);

typedef int   coor_t;    /* coordinates */
typedef char  frame_t;   /* reading frame, inc5, inc3 */
typedef float score_t;   /* scores */
typedef char  strand_t;  /* strand */

extern const coor_t   UNDEFINED_COOR;
extern const frame_t  UNDEFINED_FRAME;
extern const score_t  MIN_SCORE;
extern const score_t  MAX_SCORE;
extern const strand_t UNDEFINED_STRAND;

void     zoeCoor2Text   (coor_t, char *);
coor_t   zoeText2Coor   (const char *);
void     zoeFrame2Text  (frame_t, char *);
frame_t  zoeText2Frame  (const char *);
void     zoeScore2Text  (score_t, char *);
score_t  zoeText2Score  (const char *);
void     zoeStrand2Text (strand_t, char *);
strand_t zoeText2Strand (const char *);

void zoeS (FILE *, const char *, ...);
void zoeO (const char *, ...);
void zoeE (const char *, ...);
void zoeOs (int, ...);
void zoeEs (int, ...);
void zoeM (FILE *, int, ...);
void zoeWarn (const char *, ...);
void zoeExit  (const char *, ...);

void * zoeMalloc (size_t);
void * zoeCalloc (size_t, size_t);
void * zoeRealloc (void *, size_t);
void   zoeFree (void *);

int zoeIcmp(const void *, const void *);
int zoeFcmp(const void *, const void *);
int zoeTcmp(const void *, const void *);

struct zoeIVec  {
	int * elem;
	int   size;
	int   limit;
	int   last;
};
typedef struct zoeIVec * zoeIVec;
void    zoeDeleteIVec (zoeIVec);
zoeIVec zoeNewIVec (void);
void    zoePushIVec (zoeIVec, int);

struct zoeFVec  {
	float * elem;
	int     size;
	int     limit;
	float   last;
};
typedef struct zoeFVec * zoeFVec;
void    zoeDeleteFVec (zoeFVec);
zoeFVec zoeNewFVec (void);
void    zoePushFVec (zoeFVec, float);

struct zoeTVec  {
	char ** elem;
	int     size;
	int     limit;
	char  * last;
};
typedef struct zoeTVec * zoeTVec;
void    zoeDeleteTVec (zoeTVec);
zoeTVec zoeNewTVec (void);
void    zoePushTVec (zoeTVec, const char *);

struct zoeVec  {
	void ** elem;
	int     size;
	int     limit;
	void  * last;
};
typedef struct zoeVec * zoeVec;
void   zoeDeleteVec (zoeVec);
zoeVec zoeNewVec (void);
void   zoePushVec (zoeVec, void *);

struct zoeHash  {
	int      level;
	int      slots;
	zoeTVec  keys;
	zoeVec   vals;
	zoeVec * key;
	zoeVec * val;
};
typedef struct zoeHash * zoeHash;
void    zoeDeleteHash (zoeHash);
zoeHash zoeNewHash (void);
void    zoeSetHash (zoeHash, const char *, void *);
void *  zoeGetHash (const zoeHash, const char *);
zoeTVec zoeKeysOfHash (const zoeHash);
zoeVec  zoeValsOfHash (const zoeHash);
void    zoeStatHash (const zoeHash);

struct zoeXnode {
	zoeVec   children;
	void   * data;
	char     c;
	
};
typedef struct zoeXnode * zoeXnode;
void     zoeDeleteXnode (zoeXnode);
zoeXnode zoeNewXnode (char);
zoeXnode zoeSearchXnode (const zoeXnode, char);

struct zoeXtree {
	zoeXnode head[256];
	zoeVec   alloc;
};
typedef struct zoeXtree * zoeXtree;
void     zoeDeleteXtree (zoeXtree);
zoeXtree zoeNewXtree (void);
void *   zoeGetXtree (const zoeXtree, const char *);
void     zoeSetXtree (zoeXtree, const char *, void *);
void     zoeXtreeInfo (const zoeXtree);

struct zoeFile {
	int    type;
	FILE * stream;
	char   name[1024];
};
typedef struct zoeFile zoeFile;
void    zoeCloseFile (zoeFile);
zoeFile zoeOpenFile (const char *);

#endif
