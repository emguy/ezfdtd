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
 *     File Name : excitation.h
 * Last Modified : Thu 11 Oct 2012 12:10:26 AM EDT
 */

#ifndef EXCITATION_H
#define EXCITATION_H

typedef struct
{
    int total_points;
    int **source_points;
    int polarization;
    double *source_distribution;
    double *source_signal;
} InputPort;

int setup_point_sources (char* file_name);
int setup_input_ports (char* file_name);
int setup_hards (char* file_name);
int excite (int time_index);
int apply_ehards (int time_index);
int apply_hhards (int time_index);

#endif
