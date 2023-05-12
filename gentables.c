//Based on "How to use Fixed Point (16.16) Math" by Night Stalker,
//http://netwinder.osuosl.org/pub/netwinder/docs/nw/fix1FAQ.html
/* This program generates the trigonometric lookup tables.
It creates a header file, fixtables.h, which is used by fixmath.c

	The 'FIX_RAD_SCALE' value of 162.974 I used was derived from being able to fit an entire
360ø sine and cosine inside an array of 1024 fixed point numbers.  These
tables eat up 8K of RAM, so I figured 1024 was plenty.  Remember, fixed
point numbers are 4 bytes a piece:

   (1024 entries * 4 bytes per entry) = 4096 bytes * 2 tables = 8192 bytes.
*/
#include <stdio.h>
#include <math.h>

#include "../umath/fixmath.h"

int main(void)
{
	FILE *f;
	int i, pos;
	fix z;
	
	f = fopen("fixtables.h", "wt");
	
	pos = 0;
	fputs("#pragma once\n", f);
	fputs("fix SinTable[FIX_RAD_MAX] = {\n", f);
	for (i = 0; i < FIX_RAD_MAX-1; i++)
	{
		z = DOUBLE_TO_FIX(sin((double)i / FIX_RAD_SCALE));
		fprintf(f, "%10luU, ", (long unsigned int)z);
		pos += 12;
		if (pos > 70)
		{
			fputs("\n", f);
			pos = 0;
		}
	}
	z = DOUBLE_TO_FIX(sin((FIX_RAD_MAX-1) / FIX_RAD_SCALE));
	fprintf(f, "%10luU };\n\n", (long)z);

	pos = 0;
	fputs("fix CosTable[FIX_RAD_MAX] = {\n", f);
	for (i = 0; i < FIX_RAD_MAX-1; i++)
	{
		z = DOUBLE_TO_FIX(cos((double)i / FIX_RAD_SCALE));
		fprintf(f, "%10luU, ", (long unsigned int)z);
		pos += 12;
		if (pos > 70)
		{
			fputs("\n", f);
			pos = 0;
		}
	}
	z = DOUBLE_TO_FIX(cos((FIX_RAD_MAX-1) / FIX_RAD_SCALE));
	fprintf(f, "%10luU };\n\n", (long unsigned int)z);
	
	fclose(f);
	
	return 0;
}
