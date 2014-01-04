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
 *     File Name :   
 * Last Modified :
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "ezfdtd.h"


#ifndef ADD_COLOR
#define KNRM  ""
#define KRED  ""
#define KGRN  ""
#define KYEL  ""
#define KBLU  ""
#define KMAG  ""
#define KCYN  ""
#define KWHT  ""
#endif

#ifdef ADD_COLOR
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#endif

#define assert(f,...) if (!(f)) notice(__FILE__,__LINE__,__func__,__VA_ARGS__)
#define inspect(f,...) if (!(f)) {notice(__FILE__,__LINE__,__func__,__VA_ARGS__); return(0);}
#define check(f,...) if (!(f)) {notice(__FILE__,__LINE__,__func__,__VA_ARGS__); exit(EXIT_FAILURE);}

void notice(const char *file, int line, const char *func, const char *fmt, ...);
double power(double x, int n);

#endif /* TOOLS_H */
