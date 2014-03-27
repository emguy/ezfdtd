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
 *     File Name : ade.c
 * Last Modified : Wed 10 Oct 2012 07:10:46 PM EDT
 */

#include "tools.h"
#include "mem.h"
#include "domain.h"
#include "ade.h"
#include "h5io.h"

int total_poles;
Pole epoles[POLE_MAX];
int x, y, z;

int setup_ade (char* file_name)
{
    int total_lorentz;
    int total_drude;
    int total_debye;
    int total_sigma;
    int index_0;

    int pole_index;
    int i; 

    char dset_name[30];
    double ***tmp_a;
    double ***tmp_b;
    double ***tmp_c;

    int status;

    status = h5_get_attr(file_name, "materials", "number_of_debye_poles", &total_debye);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "materials", "number_of_drude_poles", &total_drude);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "materials", "number_of_lorentz_poles", &total_lorentz);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "materials", "ade_sigma", &total_sigma);
    inspect(status, "fail to get h5 attributes");

    total_poles = total_sigma + total_drude + total_lorentz + total_debye;
    if (total_poles == 0) return 1;

    /* sigma term */
    if (total_sigma == 1)
    {
        i = 0;
        epoles[i].px  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px, "fail to allocate memory for field");
        epoles[i].py  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py, "fail to allocate memory for field");
        epoles[i].pz  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz, "fail to allocate memory for field");
        epoles[i].px1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px1, "fail to allocate memory for field");
        epoles[i].py1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py1, "fail to allocate memory for field");
        epoles[i].pz1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz1, "fail to allocate memory for field");
        epoles[i].px2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px2, "fail to allocate memory for field");
        epoles[i].py2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py2, "fail to allocate memory for field");
        epoles[i].pz2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz2, "fail to allocate memory for field");
        epoles[i].c1  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c1, "fail to allocate memory for field");
        epoles[i].c2  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c2, "fail to allocate memory for field");
        epoles[i].c3  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c3, "fail to allocate memory for field");

        sprintf(dset_name, "/materials/sigma");
        tmp_a = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_a, "fail to load h5 dset");
        for (z = 0; z < total_z; z++)
            for (y = 0; y < total_y; y++)
                for (x = 0; x < total_x; x++)
                {
                    epoles[0].c2[x][y][z] = 1.0;
                    epoles[0].c3[x][y][z] = 2 * d_t * tmp_a[x][y][z];
                }
        free_mem3((void ***)tmp_a, total_x, total_y);
    }

    index_0 = total_sigma;
    /* lorentz terms */
    for (pole_index = 1; pole_index <= total_lorentz; pole_index++)
    {
        i = index_0 + pole_index - 1;
        epoles[i].px  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px, "fail to allocate memory for field");
        epoles[i].py  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py, "fail to allocate memory for field");
        epoles[i].pz  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz, "fail to allocate memory for field");
        epoles[i].px1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px1, "fail to allocate memory for field");
        epoles[i].py1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py1, "fail to allocate memory for field");
        epoles[i].pz1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz1, "fail to allocate memory for field");
        epoles[i].px2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px2, "fail to allocate memory for field");
        epoles[i].py2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py2, "fail to allocate memory for field");
        epoles[i].pz2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz2, "fail to allocate memory for field");
        epoles[i].c1  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c1, "fail to allocate memory for field");
        epoles[i].c2  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c2, "fail to allocate memory for field");
        epoles[i].c3  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c3, "fail to allocate memory for field");

        sprintf(dset_name, "/materials/lorentz_%02d_a", pole_index);
        tmp_a = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_a, "fail to load h5 dset");
        sprintf(dset_name, "/materials/lorentz_%02d_b", pole_index);
        tmp_b = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_b, "fail to load h5 dset");
        sprintf(dset_name, "/materials/lorentz_%02d_c", pole_index);
        tmp_c = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_c, "fail to load h5 dset");

        for (z = 0; z < total_z; z++)
            for (y = 0; y < total_y; y++)
                for (x = 0; x < total_x; x++)
                {
                    epoles[i].c1[x][y][z] = (4.0 - 2.0 * d_t * d_t * tmp_b[x][y][z]) / (2.0 + d_t * tmp_c[x][y][z]);
                    epoles[i].c2[x][y][z] = (-2.0 + d_t * tmp_c[x][y][z]) / (2.0 + d_t * tmp_c[x][y][z]);
                    epoles[i].c3[x][y][z] = (2.0 * d_t * d_t * tmp_a[x][y][z] * EPSILON0) / (2.0 + d_t * tmp_c[x][y][z]);
                    //printf("lorentz %02d [%03d %03d %03d]: %e, %e, %e\n", i, x, y, z, epoles[i].c1[x][y][z],epoles[i].c2[x][y][z],epoles[i].c3[x][y][z]);
                    //printf("lorentz %02d [%03d %03d %03d]: %e, %e, %e\n", i, x, y, z, tmp_a[x][y][z], tmp_b[x][y][z], tmp_c[x][y][z]);
                }

        free_mem3((void ***)tmp_a, total_x, total_y);
        free_mem3((void ***)tmp_b, total_x, total_y);
        free_mem3((void ***)tmp_c, total_x, total_y);
    }

    /* drude terms */
    index_0 = total_sigma + total_lorentz;
    for (pole_index = 1; pole_index <= total_drude; pole_index++)
    {
        i = index_0 + pole_index - 1;
        epoles[i].px  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px, "fail to allocate memory for field");
        epoles[i].py  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py, "fail to allocate memory for field");
        epoles[i].pz  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz, "fail to allocate memory for field");
        epoles[i].px1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px1, "fail to allocate memory for field");
        epoles[i].py1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py1, "fail to allocate memory for field");
        epoles[i].pz1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz1, "fail to allocate memory for field");
        epoles[i].px2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px2, "fail to allocate memory for field");
        epoles[i].py2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py2, "fail to allocate memory for field");
        epoles[i].pz2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz2, "fail to allocate memory for field");
        epoles[i].c1  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c1, "fail to allocate memory for field");
        epoles[i].c2  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c2, "fail to allocate memory for field");
        epoles[i].c3  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c3, "fail to allocate memory for field");

        sprintf(dset_name, "/materials/drude_%02d_a", pole_index);
        tmp_a = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_a, "fail to load h5 dset");
        sprintf(dset_name, "/materials/drude_%02d_c", pole_index);
        tmp_c = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_c, "fail to load h5 dset");

        for (z = 0; z < total_z; z++)
            for (y = 0; y < total_y; y++)
                for (x = 0; x < total_x; x++)
                {
                    epoles[i].c1[x][y][z] = 4.0 / (2.0 + d_t * tmp_c[x][y][z]);
                    epoles[i].c2[x][y][z] = (-2.0 + d_t * tmp_c[x][y][z]) / (2.0 + d_t * tmp_c[x][y][z]);
                    epoles[i].c3[x][y][z] = ( 2.0 * d_t * d_t * tmp_a[x][y][z] * EPSILON0) / (2.0 + d_t * tmp_c[x][y][z]);
                }
        free_mem3((void ***)tmp_a, total_x, total_y);
        free_mem3((void ***)tmp_c, total_x, total_y);
    }

    /* debye terms */
    index_0 = total_sigma + total_lorentz + total_drude;
    for (pole_index = 1; pole_index <= total_debye; pole_index++)
    {
        i = index_0 + pole_index - 1;
        epoles[i].px  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px, "fail to allocate memory for field");
        epoles[i].py  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py, "fail to allocate memory for field");
        epoles[i].pz  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz, "fail to allocate memory for field");
        epoles[i].px1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px1, "fail to allocate memory for field");
        epoles[i].py1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py1, "fail to allocate memory for field");
        epoles[i].pz1 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz1, "fail to allocate memory for field");
        epoles[i].px2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].px2, "fail to allocate memory for field");
        epoles[i].py2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].py2, "fail to allocate memory for field");
        epoles[i].pz2 = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].pz2, "fail to allocate memory for field");
        epoles[i].c1  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c1, "fail to allocate memory for field");
        epoles[i].c2  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c2, "fail to allocate memory for field");
        epoles[i].c3  = (double ***)mem3(type_double, total_x, total_y, total_z);
        inspect(epoles[i].c3, "fail to allocate memory for field");

        sprintf(dset_name, "/materials/debye_%02d_a", pole_index);
        tmp_a = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_a, "fail to load h5 dset");
        sprintf(dset_name, "/materials/debye_%02d_b", pole_index);
        tmp_b = (double ***)h5_load3(file_name, dset_name, total_x, total_y, total_z);
        inspect(tmp_b, "fail to load h5 dset");

        for (z = 0; z < total_z; z++)
            for (y = 0; y < total_y; y++)
                for (x = 0; x < total_x; x++)
                {
                    epoles[i].c1[x][y][z] = -2.0 * d_t* tmp_b[x][y][z];
                    epoles[i].c2[x][y][z] =  1.0;
                    epoles[i].c3[x][y][z] =  2.0 * d_t * tmp_a[x][y][z] * EPSILON0;
                }
        free_mem3((void ***)tmp_a, total_x, total_y);
        free_mem3((void ***)tmp_b, total_x, total_y);
    }

    return 1;
}

int is_dispersive (int x, int y, int z)
{
    int pole_index;
    for (pole_index = 0; pole_index < total_poles; pole_index++)
        if (epoles[pole_index].c3[x][y][z] > TOLERANCE)
            return 1;
    return 0;
}

double get_ade_ex (int x, int y, int z, double value_d, double value_e)
{
    int pole_index;
    for (pole_index = 0; pole_index < total_poles; pole_index++)
    {
        epoles[pole_index].px[x][y][z] 
            = epoles[pole_index].c1[x][y][z] * epoles[pole_index].px1[x][y][z] 
            + epoles[pole_index].c2[x][y][z] * epoles[pole_index].px2[x][y][z] 
            + epoles[pole_index].c3[x][y][z] * value_e;
        value_d = value_d - epoles[pole_index].px[x][y][z];
        epoles[pole_index].px2[x][y][z] = epoles[pole_index].px1[x][y][z];
        epoles[pole_index].px1[x][y][z] = epoles[pole_index].px[x][y][z];
    }                                 
    return value_d / epsilon[x][y][z];
}

double get_ade_ey (int x, int y, int z, double value_d, double value_e)
{
    int pole_index;
    for (pole_index = 0; pole_index < total_poles; pole_index++)
    {
        epoles[pole_index].py[x][y][z] 
            = epoles[pole_index].c1[x][y][z] * epoles[pole_index].py1[x][y][z] 
            + epoles[pole_index].c2[x][y][z] * epoles[pole_index].py2[x][y][z] 
            + epoles[pole_index].c3[x][y][z] * value_e;
        value_d = value_d  - epoles[pole_index].py[x][y][z];
        epoles[pole_index].py2[x][y][z] = epoles[pole_index].py1[x][y][z];
        epoles[pole_index].py1[x][y][z] = epoles[pole_index].py[x][y][z];
    }                                 
    return value_d / epsilon[x][y][z];
}

double get_ade_ez (int x, int y, int z, double value_d, double value_e)
{
    int pole_index;

    for (pole_index = 0; pole_index < total_poles; pole_index++)
    {
        epoles[pole_index].pz[x][y][z] 
            = epoles[pole_index].c1[x][y][z] * epoles[pole_index].pz1[x][y][z]
            + epoles[pole_index].c2[x][y][z] * epoles[pole_index].pz2[x][y][z]
            + epoles[pole_index].c3[x][y][z] * value_e;
        value_d = value_d - epoles[pole_index].pz[x][y][z];
        epoles[pole_index].pz2[x][y][z] = epoles[pole_index].pz1[x][y][z];
        epoles[pole_index].pz1[x][y][z] = epoles[pole_index].pz[x][y][z];
    }                                 
    return value_d / epsilon[x][y][z];
}

