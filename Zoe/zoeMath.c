/******************************************************************************\
 zoeMath.c - part of the ZOE library for genomic analysis
 
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

#ifndef ZOE_MATH_C
#define ZOE_MATH_C

#include "zoeMath.h"

const int zoePOWER[6][8] = {
	{1, 0,  0,   0,   0,    0,    0,      0},
	{1, 1,  1,   1,   1,    1,    1,      1},
	{1, 2,  4,   8,  16,   32,   64,    128},
	{1, 3,  9,  27,  81,  243,  729,   2187},
	{1, 4, 16,  64, 256, 1024, 4096,  16384},
	{1, 5, 25, 125, 625, 3125, 15625, 78125},
};

static double zoePI   = 3.1415926536;
static double zoeLOG2 = 0.6931471806;

double zoeLog2 (double x) {
	if (x <= 0) zoeExit("attempt to take log of %g", x);
	return (log(x) / zoeLOG2);
}

score_t zoeFloat2Score (double x) {
	if (x <= 0) return MIN_SCORE;
	return (score_t)(zoeLog2(x));
}
double  zoeScore2Float (score_t s)  {return pow((double)2, ((double)s));}

double zoeLnFactorial (int n) {
	double f;
	
	if (n == 0) return 0;
	
	f = (0.5 * log(2 * zoePI))
		+ ((n + 0.5) * log((double)n))
		- n
		+ 1 / (12 * n)
		- 1 / (360 * pow((double)n, (double)3));
	return f;
}

/* some distributions (add normal, linear, etc. later) */

score_t zoeScoreGeometric (double m, double x) {
	double f;
	double p = 1/ m;
	
	f = (x - 1) * log(1 - p) + log(p);
	f /= log((double)2);
	return (score_t) f;
}

score_t zoeScorePoisson (double m, double x) {
	double f;
	
	f = (x * log(m)) - m - zoeLnFactorial((int)x);
	f /= log((double)2);
	return (score_t) f;
}

/* base conversion - positive numbers only */

void zoeDecToBase (int n, int base, char * s) {
	int  i = 0, j;
	char c;
	
	if (base > 9) zoeExit("zoeDecToBase limit 9 (%d)\n", base);
	if (n < 0) zoeExit("zoeDecToBase positive only (%d)\n", n);
	
	do {
		s[i++] = n % base + '0';
	} while ((n /= base) > 0);
	s[i] = '\0';
	
	for (i = 0, j = strlen(s) -1; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

int zoeBaseToDec (int base, const char * s) {
	char c;
	int  i, len, val = 0;
	
	if (base > 9) zoeExit("zoeBaseToDec limit 9 (%d)\n", base);
	
	len = strlen(s);
	for (i = 0; i < len; i++) {
		c = s[i] - '0';
		if (c < 0 || c > 9) zoeExit("zoeBaseToDec illegal symbol (%c)", s[i]);
		val += pow((float)base, (float)(len -i -1)) * c;
	}
	
	return val;
}

double zoeDivide (double a, double b) {
	if (a && b) {
		return a / b;
	} else if (a) {
		return FLT_MAX;
	} else if (b) {
		return 0;
	} else {
		return 1;
	}
}


#endif
