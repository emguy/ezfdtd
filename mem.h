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
 *     File Name : matrix.h
 * Last Modified : Tue 24 Sep 2013 05:59:03 PM EDT
 */


#ifndef MEMORY_H
#define MEMORY_H

enum {type_null, type_int, type_double, type_float, type_char};

void ***mem3(int type, size_t xmax, size_t ymax, size_t zmax);
void **mem2(int type, size_t xmax, size_t ymax);
void *mem1(int type, size_t x);
void **contiguous_mem2 (int type, size_t xmax, size_t ymax);
void ***contiguous_mem3 (int type, size_t xmax, size_t ymax, size_t zmax);
void free_mem3(void ***mem, int xmax, int ymax);
void free_mem2(void **mem, int xmax);
void free_mem1(void *mem);
void free_contiguous_mem2 (void **mem);
void free_contiguous_mem3 (void ***mem);

#endif
