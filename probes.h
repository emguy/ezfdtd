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
 *     File Name : probes.h
 * Last Modified : Fri 12 Oct 2012 03:20:37 PM EDT
 */

#ifndef PROBES_H
#define PROBES_H

typedef struct
{
    int write_timedomain;
    int **probe_points;
    int total_points;
    int polarization;
    double **fields;
    double **dft_real;
    double **dft_imag;
} OutputPort;

typedef struct
{
    int dim;
    int slice_index;
    int polarization;
    double **field_buffer;
    int with_pml;
} FieldPlane;

extern OutputPort output_ports[100];
extern int total_output_ports;

int setup_planes(char* file_name, char* output_file_name, int with_pml);
int update_planes (char* file_name, int time_index);
int setup_output_ports(char* file_name);
int update_ports(int time_index);
int write_ports(char* file_name);

#endif
