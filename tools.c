/* NOTE: This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either 
 * version 3, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * Bugs can be reported to Yu Zhang <zhang235@mcmaster.ca>.
 *
 *     File Name : tools.c
 * Last Modified : Wed 25 Sep 2013 12:26:30 PM EDT
 */

#include "tools.h"
#define DEBUG 1

void notice(const char *file, int line, const char *func, const char *fmt, ...)
{
	va_list ap;
	va_start(ap,fmt);
	if (DEBUG)
	{
		if (file) fprintf(stderr, KRED "%s ", file);
		if (line) fprintf(stderr, KRED "%d ", line);
		if (func) fprintf(stderr, KRED "in %s() ", func);
		fprintf(stderr, KRED ": " KCYN);
	}
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fprintf(stderr, KNRM "\n");
}
