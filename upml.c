/* NOTE: This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either 
 * version 3, or (at your option) any later version.  *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * Bugs can be reported to Yu Zhang <zhang235@mcmaster.ca>.
 *
 *     File Name : upml.c
 * Last Modified : Fri 12 Oct 2012 08:21:29 PM EDT
 */
#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "domain.h"
#include "ezfdtd.h"
#include "h5in.h"
#include "upml.h"

double **pml_c1;
double **pmlC2;
double *boundaryEpsilon;

void setup_pml (char file_name[])
{
    double reflectionCoefficient0; 
    int grading_order; 

    double dh;
    double boundaryThickness;
    double pml_sigma_hax;
    double pml_factor;
    double pmlSigma;

    int layer_index, partition_index;
    char attr_name[30];

    reflectionCoefficient0 = getDoubleAttribute(file_name, "boundaries", "PML_reflectionCoefficient0");
    grading_order = getDoubleAttribute(file_name, "boundaries", "PML_grading_order");

    boundaryEpsilon = allocateMemory1d(7, 1.0);
    pml_c1 = allocateMemory(abc_size, 7, 1.0);
    pmlC2 = allocateMemory(abc_size, 7, d_t);

    for (partition_index = 1; partition_index < 7; partition_index++)
    {
        switch (partition_index)
        {
            case 1: dh = d_z; break;
            case 2: dh = d_z; break;
            case 3: dh = d_x; break;
            case 4: dh = d_x; break;
            case 5: dh = d_y; break;
            case 6: dh = d_y; break;
            default: dh = d_z;
        }
        sprintf(attr_name, "boundary_epsilon_%01d", partition_index);
        boundaryEpsilon[partition_index] = getDoubleAttribute(file_name, "boundaries", "boundary_1_epsilon");

        boundaryThickness = (double )abc_size * dh;
        pml_sigma_hax = -log(reflectionCoefficient0) * (grading_order + 1.0) 
            * boundaryEpsilon[partition_index] * C0 / (2.0 * boundaryThickness);
        pml_factor = pml_sigma_hax 
            / ((grading_order + 1) * pow(2, grading_order + 1) * pow(abc_size, grading_order));

        for (layer_index = 0; layer_index < abc_size; layer_index++)
        {
            if (layer_index == 0)
                pmlSigma = pml_factor;
            else
                pmlSigma = pml_factor * (pow(2 * layer_index + 1, grading_order + 1) 
                        - pow(2 * layer_index - 1, grading_order + 1));

            pml_c1[layer_index][partition_index] = (2 * boundaryEpsilon[partition_index] * EPSILON0 - pmlSigma * d_t)
                / (2 * boundaryEpsilon[partition_index] * EPSILON0 + pmlSigma * d_t);
            pmlC2[layer_index][partition_index] = (2 * boundaryEpsilon[partition_index] * EPSILON0 * d_t)
                / (2 * boundaryEpsilon[partition_index] * EPSILON0 + pmlSigma * d_t);
            if (partition_index == 5)
                printf("%e, %e\n", pml_c1[layer_index][partition_index], pmlC2[layer_index][partition_index]);
        }
    }
}

void evaluatePmlH (int x, int y, int z)
{
    int partition_index;
    int layer_index;
    double tmp;

    partition_index = get_partition(x, y, z);

    if (partition_index != 0 && partition_data[partition_index].boundary_type == 3)
    {
        switch (partition_index)
        {
            case 1: layer_index = partition_data[partition_z0].z_stop - z - 1; break;
            case 2: layer_index = z -partition_data[partition_z1].z_start;     break;
            case 3: layer_index = partition_data[partition_x0].x_stop - x - 1; break;
            case 4: layer_index = x - partition_data[partition_x1].x_start;    break;
            case 5: layer_index = partition_data[partition_y0].y_stop - y - 1; break;
            case 6: layer_index = y - partition_data[partition_y1].y_start;    break;
            default: layer_index = 0;
        }

        tmp = bx[x][y][z];
        bx[x][y][z] = pml_c1[layer_index][partition_index] * bx[x][y][z] 
            -pmlC2[layer_index][partition_index] * (ez[x][y+1][z] - ez[x][y][z]) / d_y
            +pmlC2[layer_index][partition_index] * (ey[x][y][z+1] - ey[x][y][z]) / d_z;
        hx[x][y][z] = pml_c1[layer_index][partition_index] * hx[x][y][z] + bx[x][y][z] / MU0
            - pml_c1[layer_index][partition_index] * tmp / MU0;

        tmp = by[x][y][z];
        by[x][y][z] = pml_c1[layer_index][partition_index] * by[x][y][z] 
            -pmlC2[layer_index][partition_index] * (ex[x][y+1][z] - ex[x][y][z]) / d_z
            +pmlC2[layer_index][partition_index] * (ez[x][y][z+1] - ez[x][y][z]) / d_x;
        hy[x][y][z] = pml_c1[layer_index][partition_index] * hy[x][y][z] + by[x][y][z] / MU0
            - pml_c1[layer_index][partition_index] * tmp / MU0;

        tmp = bz[x][y][z];
        bz[x][y][z] = pml_c1[layer_index][partition_index] * bz[x][y][z] 
            -pmlC2[layer_index][partition_index] * (ey[x][y+1][z] - ey[x][y][z]) / d_x
            +pmlC2[layer_index][partition_index] * (ex[x][y][z+1] - ex[x][y][z]) / d_y;
        hz[x][y][z] = pml_c1[layer_index][partition_index] * hz[x][y][z] + bz[x][y][z] / MU0
            - pml_c1[layer_index][partition_index] * tmp / MU0;
    }
}

void evaluatePmlE (int x, int y, int z)
{
    int partition_index;
    int layer_index;
    double tmp;

    partition_index = get_partition(x, y, z);

    if (partition_index != 0 && partition_data[partition_index].boundary_type == 3)
    {
        switch (partition_index)
        {
            case 1: layer_index = partition_data[partition_z0].z_stop - z - 1; break;
            case 2: layer_index =-partition_data[partition_z1].z_start + z;    break;
            case 3: layer_index = partition_data[partition_x0].x_stop - x - 1; break;
            case 4: layer_index =-partition_data[partition_x1].x_start + x;    break;
            case 5: layer_index = partition_data[partition_y0].y_stop - y - 1; break;
            case 6: layer_index =-partition_data[partition_y1].y_start + y;    break;
            default: layer_index = 0;
        }
        tmp = dx[x][y][z];
        dx[x][y][z] = pml_c1[layer_index][partition_index] * dx[x][y][z] 
            +pmlC2[layer_index][partition_index] * (hz[x][y+1][z] - hz[x][y][z]) / d_y
            -pmlC2[layer_index][partition_index] * (hy[x][y][z+1] - hy[x][y][z]) / d_z;
        ex[x][y][z] = pml_c1[layer_index][partition_index] * ex[x][y][z] 
            + dx[x][y][z] / boundaryEpsilon[partition_index]
            - pml_c1[layer_index][partition_index] * tmp / boundaryEpsilon[partition_index];

        tmp = dy[x][y][z];
        dy[x][y][z] = pml_c1[layer_index][partition_index] * dy[x][y][z] 
            +pmlC2[layer_index][partition_index] * (hx[x][y+1][z] - hx[x][y][z]) / d_z
            -pmlC2[layer_index][partition_index] * (hz[x][y][z+1] - hz[x][y][z]) / d_x;
        ey[x][y][z] = pml_c1[layer_index][partition_index] * ey[x][y][z] 
            + dy[x][y][z] / boundaryEpsilon[partition_index]
            - pml_c1[layer_index][partition_index] * tmp / boundaryEpsilon[partition_index];

        tmp = dz[x][y][z];
        dz[x][y][z] = pml_c1[layer_index][partition_index] * dz[x][y][z] 
            +pmlC2[layer_index][partition_index] * (hy[x][y+1][z] - hy[x][y][z]) / d_x
            -pmlC2[layer_index][partition_index] * (hx[x][y][z+1] - hx[x][y][z]) / d_y;
        ez[x][y][z] = pml_c1[layer_index][partition_index] * ez[x][y][z] 
            + dz[x][y][z] / boundaryEpsilon[partition_index]
            - pml_c1[layer_index][partition_index] * tmp / boundaryEpsilon[partition_index];
    }
}
