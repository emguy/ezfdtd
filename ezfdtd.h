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
 *     File Name : ezfdtd.h
 * Last Modified : Sat 06 Oct 2012 11:16:54 PM EDT
 */

#ifndef EZFDTD_H
#define EZFDTD_H

#ifndef DEBUG
#define DEBUG 1
#endif

#define POLE_MAX 10
#define TOLERANCE 6.6e-66
#define PI 3.1415926535897932384626434


typedef struct _DomainData
{
    unsigned int x_start; 
    unsigned int y_start;
    unsigned int z_start;
    unsigned int x_stop;
    unsigned int y_stop;
    unsigned int z_stop;
    unsigned int boundary_type;
    unsigned int thickness;
    unsigned int size_x;
    unsigned int size_y;
    unsigned int size_z;
} DomainData;

/* dimension */
enum{dim_t, dim_x, dim_y, dim_z};

/* partition index */
enum{partition_main, partition_z0, partition_z1, partition_x0, partition_x1, partition_y0, partition_y1};

/* boundary types */
enum{boundary_air, boundary_pec, boundary_pmc, boundary_pml, boundary_mur};

/* mode */
extern unsigned int mode;
enum{mode_full, mode_tmx, mode_tmy, mode_tmz, mode_tex, mode_tey, mode_tez};

/* polarization */
enum{p_0, p_ex, p_ey, p_ez, p_hx, p_hy, p_hz};

/* pml types */
extern unsigned int pml_type;
enum{pml_pml, pml_cpml};

/* constants */
extern double MU0;
extern double C0;
extern double EPSILON0;

/* discretization */
extern double d_x;
extern double d_y;
extern double d_z;
extern double d_t;

/* simulation time */
extern unsigned int total_timesteps;

/* grid size */
extern unsigned int abc_size;
extern unsigned int total_length_x;
extern unsigned int total_length_y; 
extern unsigned int total_length_z; 
extern unsigned int main_length_x;
extern unsigned int main_length_y;
extern unsigned int main_length_z;

/* the fields */
extern double ***ex;
extern double ***ey;
extern double ***ez;
extern double ***hx;
extern double ***hy;
extern double ***hz;
extern double ***dx;
extern double ***dy;
extern double ***dz;
extern double ***bx;
extern double ***by;
extern double ***bz;

/* the excitation sources */
extern double ***dipole_ex;
extern double ***dipole_ey;
extern double ***dipole_ez;
extern double ***dipole_hx;
extern double ***dipole_hy;
extern double ***dipole_hz;

/* domain partition */
extern DomainData partition_data[7];

#endif
