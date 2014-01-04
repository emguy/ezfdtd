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
 *     File Name : ade.h
 * Last Modified : Wed 10 Oct 2012 07:28:19 PM EDT
 */

#ifndef ADE_H
#define ADE_H

typedef struct
{
    int type;
    int ***region;
    double ***c1;
    double ***c2;
    double ***c3;
    double ***px;
    double ***py;
    double ***pz;
    double ***px1;
    double ***py1;
    double ***pz1;
    double ***px2;
    double ***py2;
    double ***pz2;
}Pole;

int setup_ade (char* file_name);
int is_dispersive (int x_main, int y_main, int z_main);
double get_ade_ex (int x_main, int y_main, int z_main, double value_d, double value_e);
double get_ade_ey (int x_main, int y_main, int z_main, double value_d, double value_e);
double get_ade_ez (int x_main, int y_main, int z_main, double value_d, double value_e);

#endif
