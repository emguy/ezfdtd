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
 *     File Name : cpml.h
 * Last Modified : Mon 03 Dec 2012 01:33:10 PM EST
 */


#ifndef CPML_H
#define CPML_H
typedef struct
{
    double ***exy;
    double ***exz;
    double ***eyx;
    double ***eyz;
    double ***ezx;
    double ***ezy;
    double ***hxy;
    double ***hxz;
    double ***hyx;
    double ***hyz;
    double ***hzx;
    double ***hzy;
} CPMLFields;


int setup_cpml (char* file_name);
void cpml_get_d();
void cpml_get_b();

#endif
