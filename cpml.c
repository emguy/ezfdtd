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
 *     File Name : cpml.c
 * Last Modified : Thu 06 Dec 2012 10:39:02 AM EST
 */

#include <math.h>
#include "tools.h"
#include "mem.h"
#include "domain.h"
#include "h5io.h"
#include "cpml.h"

double *pml_bex;
double *pml_bey;
double *pml_bez;
double *pml_bhx;
double *pml_bhy;
double *pml_bhz;

double *pml_cex;
double *pml_cey;
double *pml_cez;
double *pml_chx;
double *pml_chy;
double *pml_chz;

double *pml_kappa_h;
double *pml_kappa_e;

CPMLFields psi[7];

double bex;
double bey;
double bez;
double bhx;
double bhy;
double bhz;

double cex;
double cey;
double cez;
double chx;
double chy;
double chz;

double kappa_h;
double kappa_e;

int setup_cpml (char* file_name)
{
    double pml_grading_order; 
    double pml_sigma_max;
    double pml_alpha_max;
    double pml_kappa_max;
    double pml_sigma_e;
    double pml_sigma_h;
    double pml_alpha_e;
    double pml_alpha_h;

    int layer_index, partition_index;
    int x_start, x_stop, y_start, y_stop, z_start, z_stop;

    int status;

    status = h5_get_attr(file_name, "boundaries", "pml_sigma_max", &pml_sigma_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_alpha_max", &pml_alpha_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_kappa_max", &pml_kappa_max);
    inspect(status, "fail to get h5 attributes");
    status = h5_get_attr(file_name, "boundaries", "pml_grading_order", &pml_grading_order);
    inspect(status, "fail to get h5 attributes");
    //printf("CPML: %e %e %e %e\n", pml_sigma_max, pml_alpha_max, pml_kappa_max, pml_grading_order);


    pml_bex  = (double *)mem1(type_double, abc_size);
    inspect(pml_bex, "fail to allocate memory for field");
    pml_bey  = (double *)mem1(type_double, abc_size);
    inspect(pml_bey, "fail to allocate memory for field");
    pml_bez  = (double *)mem1(type_double, abc_size);
    inspect(pml_bez, "fail to allocate memory for field");
    pml_bhx  = (double *)mem1(type_double, abc_size);
    inspect(pml_bhx, "fail to allocate memory for field");
    pml_bhy  = (double *)mem1(type_double, abc_size);
    inspect(pml_bhy, "fail to allocate memory for field");
    pml_bhz  = (double *)mem1(type_double, abc_size);
    inspect(pml_bhz, "fail to allocate memory for field");
    pml_cex = (double *)mem1(type_double, abc_size);
    inspect(pml_cex, "fail to allocate memory for field");
    pml_cey = (double *)mem1(type_double, abc_size);
    inspect(pml_cey, "fail to allocate memory for field");
    pml_cez = (double *)mem1(type_double, abc_size);
    inspect(pml_cez, "fail to allocate memory for field");
    pml_chx = (double *)mem1(type_double, abc_size);
    inspect(pml_chx, "fail to allocate memory for field");
    pml_chy = (double *)mem1(type_double, abc_size);
    inspect(pml_chy, "fail to allocate memory for field");
    pml_chz = (double *)mem1(type_double, abc_size);
    inspect(pml_chz, "fail to allocate memory for field");
    pml_kappa_e = (double *)mem1(type_double, abc_size);
    inspect(pml_kappa_e, "fail to allocate memory for field");
    pml_kappa_h = (double *)mem1(type_double, abc_size);
    inspect(pml_kappa_h, "fail to allocate memory for field");

    for (layer_index = 0; layer_index < abc_size; layer_index++)
    {
        pml_sigma_e = pml_sigma_max * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_sigma_h = pml_sigma_max * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_kappa_e[layer_index] = 1.0 + (pml_kappa_max - 1.0) * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_kappa_h[layer_index] = 1.0 + (pml_kappa_max - 1.0) * pow((layer_index + 1.0) / (abc_size + 1.0), pml_grading_order);
        pml_alpha_e = pml_alpha_max * (abc_size - layer_index - 1.0) / (abc_size - 1.0);
        pml_alpha_h = pml_alpha_max * (abc_size - layer_index - 1.0) / (abc_size - 1.0);

        //printf("layer %03d: %e, %e, %e, %e \n", layer_index, pml_sigma_e, pml_sigma_h, pml_alpha_e, pml_alpha_h);
        //printf("layer %03d: %e\n", layer_index, pml_kappa_e);

        pml_bex[layer_index]  = exp(-(pml_sigma_e / pml_kappa_e[layer_index] + pml_alpha_e) * d_t / EPSILON0);
        pml_bey[layer_index]  = exp(-(pml_sigma_e / pml_kappa_e[layer_index] + pml_alpha_e) * d_t / EPSILON0);
        pml_bez[layer_index]  = exp(-(pml_sigma_e / pml_kappa_e[layer_index] + pml_alpha_e) * d_t / EPSILON0);
        pml_bhx[layer_index]  = exp(-(pml_sigma_h / pml_kappa_h[layer_index] + pml_alpha_h) * d_t / EPSILON0);
        pml_bhy[layer_index]  = exp(-(pml_sigma_h / pml_kappa_h[layer_index] + pml_alpha_h) * d_t / EPSILON0);
        pml_bhz[layer_index]  = exp(-(pml_sigma_h / pml_kappa_h[layer_index] + pml_alpha_h) * d_t / EPSILON0);

        pml_cex[layer_index] = pml_sigma_e * (pml_bex[layer_index] - 1) / (pml_sigma_e + pml_alpha_e * pml_kappa_e[layer_index]) / d_x / pml_kappa_e[layer_index];
        pml_cey[layer_index] = pml_sigma_e * (pml_bey[layer_index] - 1) / (pml_sigma_e + pml_alpha_e * pml_kappa_e[layer_index]) / d_y / pml_kappa_e[layer_index];
        pml_cez[layer_index] = pml_sigma_e * (pml_bez[layer_index] - 1) / (pml_sigma_e + pml_alpha_e * pml_kappa_e[layer_index]) / d_z / pml_kappa_e[layer_index];

        pml_chx[layer_index] = pml_sigma_h * (pml_bex[layer_index] - 1) / (pml_sigma_h + pml_alpha_h * pml_kappa_h[layer_index]) / d_x / pml_kappa_h[layer_index];
        pml_chy[layer_index] = pml_sigma_h * (pml_bey[layer_index] - 1) / (pml_sigma_h + pml_alpha_h * pml_kappa_h[layer_index]) / d_y / pml_kappa_h[layer_index];
        pml_chz[layer_index] = pml_sigma_h * (pml_bez[layer_index] - 1) / (pml_sigma_h + pml_alpha_h * pml_kappa_h[layer_index]) / d_z / pml_kappa_h[layer_index];

        //printf("layer %03d: %e, %e, %e, %e, %e, %e, %e  \n", layer_index, pml_sigma_e, pml_kappa_e[layer_index], pml_alpha_e, be[layer_index], pml_cex[layer_index], bh[layer_index], pml_chx[layer_index]);
    }

    for (partition_index = 1; partition_index < 7; partition_index++)
    {
        x_start = partition_data[partition_index].x_start;
        x_stop = partition_data[partition_index].x_stop;
        y_start = partition_data[partition_index].y_start;
        y_stop = partition_data[partition_index].y_stop;
        z_start = partition_data[partition_index].z_start;
        z_stop = partition_data[partition_index].z_stop;

        if (partition_data[partition_index].boundary_type == boundary_pml)
        {
            switch (mode)
            {
                case mode_tmx:
                    psi[partition_index].hyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyz, "fail to allocate memory for field");
                    psi[partition_index].hzy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzy, "fail to allocate memory for field");
                    psi[partition_index].eyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyx, "fail to allocate memory for field");
                    psi[partition_index].ezx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezx, "fail to allocate memory for field");
                    break;
                case mode_tmy:
                    psi[partition_index].hxz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxz, "fail to allocate memory for field");
                    psi[partition_index].hzx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzx, "fail to allocate memory for field");
                    psi[partition_index].exy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exy, "fail to allocate memory for field");
                    psi[partition_index].ezy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezy, "fail to allocate memory for field");
                    break;
                case mode_tmz:
                    psi[partition_index].hyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyx, "fail to allocate memory for field");
                    psi[partition_index].hxy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxy, "fail to allocate memory for field");
                    psi[partition_index].eyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyz, "fail to allocate memory for field");
                    psi[partition_index].exz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exz, "fail to allocate memory for field");
                    break;
                case mode_tex:
                    psi[partition_index].hyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyx, "fail to allocate memory for field");
                    psi[partition_index].hzx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzx, "fail to allocate memory for field");
                    psi[partition_index].eyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyz, "fail to allocate memory for field");
                    psi[partition_index].ezy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezy, "fail to allocate memory for field");
                    break;
                case mode_tey:
                    psi[partition_index].hxy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxy, "fail to allocate memory for field");
                    psi[partition_index].hzy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzy, "fail to allocate memory for field");
                    psi[partition_index].exz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exz, "fail to allocate memory for field");
                    psi[partition_index].ezx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezx, "fail to allocate memory for field");
                    break;
                case mode_tez:
                    psi[partition_index].hxz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxz, "fail to allocate memory for field");
                    psi[partition_index].hyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyz, "fail to allocate memory for field");
                    psi[partition_index].exy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exy, "fail to allocate memory for field");
                    psi[partition_index].eyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyx, "fail to allocate memory for field");
                    break;
                case mode_full:
                    psi[partition_index].hxy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxy, "fail to allocate memory for field");
                    psi[partition_index].hxz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hxz, "fail to allocate memory for field");
                    psi[partition_index].hyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyx, "fail to allocate memory for field");
                    psi[partition_index].hyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hyz, "fail to allocate memory for field");
                    psi[partition_index].hzx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzx, "fail to allocate memory for field");
                    psi[partition_index].hzy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].hzy, "fail to allocate memory for field");
                    psi[partition_index].exy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exy, "fail to allocate memory for field");
                    psi[partition_index].exz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].exz, "fail to allocate memory for field");
                    psi[partition_index].eyx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyx, "fail to allocate memory for field");
                    psi[partition_index].eyz = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].eyz, "fail to allocate memory for field");
                    psi[partition_index].ezx = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezx, "fail to allocate memory for field");
                    psi[partition_index].ezy = (double ***)mem3(type_double, x_stop - x_start, y_stop - y_start, z_stop - z_start);
                    inspect(psi[partition_index].ezy, "fail to allocate memory for field");
                    break;
            }
        }
    }
    return 1;
}

void cpml_get_d()
{
    int x, y, z;
    int partition_index;
    int x_start;
    int y_start; 
    int z_start;

    for (x = 0; x < total_length_x; x++)
    {
        for (y = 0; y < total_length_y; y++)
        {
            for (z = 0; z < total_length_z; z++)
            {
                if (in_partition_main(x, y, z)) continue;

                partition_index = get_partition(x, y, z);
                x_start = partition_data[partition_index].x_start;
                y_start = partition_data[partition_index].y_start;
                z_start = partition_data[partition_index].z_start;

                bex = 1.0;
                bey = 1.0;
                bez = 1.0;
                cex = 0.0;
                cey = 0.0;
                cez = 0.0;

                if (in_partition_x0(x, y, z))
                {
                    bex = pml_bex[partition_data[partition_x0].x_stop - x - 1];
                    cex = pml_cex[partition_data[partition_x0].x_stop - x - 1];
                    kappa_e = pml_kappa_e[partition_data[partition_x0].x_stop - x - 1];
                }
                if (in_partition_y0(x, y, z))
                {
                    bey = pml_bey[partition_data[partition_y0].y_stop - y - 1];
                    cey = pml_cey[partition_data[partition_y0].y_stop - y - 1];
                    kappa_e = pml_kappa_e[partition_data[partition_y0].y_stop - y - 1];
                }
                if (in_partition_z0(x, y, z))
                {
                    bez = pml_bez[partition_data[partition_z0].z_stop - z - 1];
                    cez = pml_cez[partition_data[partition_z0].z_stop - z - 1];
                    kappa_e = pml_kappa_e[partition_data[partition_z0].z_stop - z - 1];
                }
                if (in_partition_x1(x, y, z))
                {
                    bex = pml_bex[-partition_data[partition_x1].x_start + x];
                    cex = pml_cex[-partition_data[partition_x1].x_start + x];
                    kappa_e = pml_kappa_e[-partition_data[partition_x1].x_start + x];
                }
                if (in_partition_y1(x, y, z))
                {
                    bey = pml_bey[-partition_data[partition_y1].y_start + y];
                    cey = pml_cey[-partition_data[partition_y1].y_start + y];
                    kappa_e = pml_kappa_e[-partition_data[partition_y1].y_start + y];
                }
                if (in_partition_z1(x, y, z))
                {
                    bez = pml_bez[-partition_data[partition_z1].z_start + z];
                    cez = pml_cez[-partition_data[partition_z1].z_start + z];
                    kappa_e = pml_kappa_e[-partition_data[partition_z1].z_start + z];
                }

                if (psi[partition_index].hxy && x != 0)
                    psi[partition_index].hxy[x-x_start][y-y_start][z-z_start] 
                        = bex * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start]
                        + cex * (hy[x][y][z] - hy[x-1][y][z]);
                if (psi[partition_index].hxz && x != 0)
                    psi[partition_index].hxz[x-x_start][y-y_start][z-z_start] 
                        = bex * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start]
                        + cex * (hz[x][y][z] - hz[x-1][y][z]);
                if (psi[partition_index].hyx && y != 0)
                    psi[partition_index].hyx[x-x_start][y-y_start][z-z_start] 
                        = bey * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start]
                        + cey * (hx[x][y][z] - hx[x][y-1][z]);
                if (psi[partition_index].hyz && y != 0)
                    psi[partition_index].hyz[x-x_start][y-y_start][z-z_start] 
                        = bey * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start]
                        + cey * (hz[x][y][z] - hz[x][y-1][z]);
                if (psi[partition_index].hzx && z != 0)
                    psi[partition_index].hzx[x-x_start][y-y_start][z-z_start] 
                        = bez * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start]
                        + cez * (hx[x][y][z] - hx[x][y][z-1]);
                if (psi[partition_index].hzy && z != 0)
                    psi[partition_index].hzy[x-x_start][y-y_start][z-z_start] 
                        = bez * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start]
                        + cez * (hy[x][y][z] - hy[x][y][z-1]);

                switch (mode)
                {
                    case mode_full:
                        if (y != 0)
                            dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]) / kappa_e;
                        if (z != 0)
                            dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]) / kappa_e;
                        if (z != 0)
                            dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]) / kappa_e;
                        if (x != 0)
                            dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]) / kappa_e;
                        if (x != 0)
                            dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]) / kappa_e;
                        if (y != 0)
                            dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]) / kappa_e;
                        break;
                    case mode_tmx:
                        if (y != 0)
                            dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]) / kappa_e;
                        if (z != 0)
                            dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]) / kappa_e;
                        break;
                    case mode_tmy:
                        if (z != 0)
                            dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]) / kappa_e;
                        if (x != 0)
                            dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]) / kappa_e;
                        break;
                    case mode_tmz:
                        if (x != 0)
                            dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]) / kappa_e;
                        if (y != 0)
                            dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]) / kappa_e;
                        break;
                    case mode_tex:
                        if (z != 0)
                            dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]) / kappa_e;
                        if (y != 0)
                            dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]) / kappa_e;
                        break;
                    case mode_tey:
                        if (z != 0)
                            dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]) / kappa_e;
                        if (x != 0)
                            dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]) / kappa_e;
                        break;
                    case mode_tez:
                        if (y != 0)
                            dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]) / kappa_e;
                        if (x != 0)
                            dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]) / kappa_e;
                        break;
                }

                switch (mode)
                {
                    case mode_full:
                        dx[x][y][z] = dx[x][y][z] + d_t * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start];
                        dx[x][y][z] = dx[x][y][z] - d_t * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start];
                        dy[x][y][z] = dy[x][y][z] + d_t * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start];
                        dy[x][y][z] = dy[x][y][z] - d_t * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start];
                        dz[x][y][z] = dz[x][y][z] + d_t * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start];
                        dz[x][y][z] = dz[x][y][z] - d_t * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tmx:
                        dx[x][y][z] = dx[x][y][z] + d_t * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start];
                        dx[x][y][z] = dx[x][y][z] - d_t * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tmy:
                        dy[x][y][z] = dy[x][y][z] + d_t * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start];
                        dy[x][y][z] = dy[x][y][z] - d_t * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tmz:
                        dz[x][y][z] = dz[x][y][z] + d_t * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start];
                        dz[x][y][z] = dz[x][y][z] - d_t * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tex:
                        dy[x][y][z] = dy[x][y][z] + d_t * psi[partition_index].hzx[x-x_start][y-y_start][z-z_start];
                        dz[x][y][z] = dz[x][y][z] - d_t * psi[partition_index].hyx[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tey:
                        dx[x][y][z] = dx[x][y][z] - d_t * psi[partition_index].hzy[x-x_start][y-y_start][z-z_start];
                        dz[x][y][z] = dz[x][y][z] + d_t * psi[partition_index].hxy[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tez:
                        dx[x][y][z] = dx[x][y][z] + d_t * psi[partition_index].hyz[x-x_start][y-y_start][z-z_start];
                        dy[x][y][z] = dy[x][y][z] - d_t * psi[partition_index].hxz[x-x_start][y-y_start][z-z_start];
                        break;
                }
            }
        }
    }
}


void cpml_get_b ()
{
    int x, y, z;
    int partition_index;
    int x_start;
    int y_start; 
    int z_start;

    for (x = 0; x < total_length_x; x++)
    {
        for (y = 0; y < total_length_y; y++)
        {
            for (z = 0; z < total_length_z; z++)
            {
                if (in_partition_main(x, y, z)) continue;

                partition_index = get_partition(x, y, z);
                x_start = partition_data[partition_index].x_start;
                y_start = partition_data[partition_index].y_start;
                z_start = partition_data[partition_index].z_start;

                bhx = 1.0;
                bhy = 1.0;
                bhz = 1.0;
                chx = 0.0;
                chy = 0.0;
                chz = 0.0;

                if (in_partition_x0(x, y, z))
                {
                    bhx = pml_bhx[partition_data[partition_x0].x_stop - x - 1];
                    chx = pml_chx[partition_data[partition_x0].x_stop - x - 1];
                    kappa_h = pml_kappa_h[partition_data[partition_x0].x_stop - x - 1];
                }
                if (in_partition_y0(x, y, z))
                {
                    bhy = pml_bhy[partition_data[partition_y0].y_stop - y - 1];
                    chy = pml_chy[partition_data[partition_y0].y_stop - y - 1];
                    kappa_h = pml_kappa_h[partition_data[partition_y0].y_stop - y - 1];
                }
                if (in_partition_z0(x, y, z))
                {
                    bhz = pml_bhz[partition_data[partition_z0].z_stop - z - 1];
                    chz = pml_chz[partition_data[partition_z0].z_stop - z - 1];
                    kappa_h = pml_kappa_h[partition_data[partition_z0].z_stop - z - 1];
                }
                if (in_partition_x1(x, y, z))
                {
                    bhx = pml_bhx[-partition_data[partition_x1].x_start + x];
                    chx = pml_chx[-partition_data[partition_x1].x_start + x];
                    kappa_h = pml_kappa_h[-partition_data[partition_x1].x_start + x];
                }
                if (in_partition_y1(x, y, z))
                {
                    bhy = pml_bhy[-partition_data[partition_y1].y_start + y];
                    chy = pml_chy[-partition_data[partition_y1].y_start + y];
                    kappa_h = pml_kappa_h[-partition_data[partition_y1].y_start + y];
                }
                if (in_partition_z1(x, y, z))
                {
                    bhz = pml_bhz[-partition_data[partition_z1].z_start + z];
                    chz = pml_chz[-partition_data[partition_z1].z_start + z];
                    kappa_h = pml_kappa_h[-partition_data[partition_z1].z_start + z];
                }

                if (psi[partition_index].exy)
                    psi[partition_index].exy[x-x_start][y-y_start][z-z_start] 
                        = bhx * psi[partition_index].exy[x-x_start][y-y_start][z-z_start]
                        + chx * (ey[x+1][y][z] - ey[x][y][z]);
                if (psi[partition_index].exz)
                    psi[partition_index].exz[x-x_start][y-y_start][z-z_start] 
                        = bhx * psi[partition_index].exz[x-x_start][y-y_start][z-z_start]
                        + chx * (ez[x+1][y][z] - ez[x][y][z]);
                if (psi[partition_index].eyx)
                    psi[partition_index].eyx[x-x_start][y-y_start][z-z_start] 
                        = bhy * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start]
                        + chy * (ex[x][y+1][z] - ex[x][y][z]);
                if (psi[partition_index].eyz)
                    psi[partition_index].eyz[x-x_start][y-y_start][z-z_start] 
                        = bhy * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start]
                        + chy * (ez[x][y+1][z] - ez[x][y][z]);
                if (psi[partition_index].ezx)
                    psi[partition_index].ezx[x-x_start][y-y_start][z-z_start] 
                        = bhy * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start]
                        + chz * (ex[x][y][z+1] - ex[x][y][z]);
                if (psi[partition_index].ezy)
                    psi[partition_index].ezy[x-x_start][y-y_start][z-z_start] 
                        = bhy * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start]
                        + chz * (ey[x][y][z+1] - ey[x][y][z]);

                switch (mode)
                {
                    case mode_full:
                        bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]) / kappa_h;
                        bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]) / kappa_h;
                        by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]) / kappa_h;
                        by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]) / kappa_h;
                        bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]) / kappa_h;
                        bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]) / kappa_h;
                        break;
                    case mode_tmx:
                        by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]) / kappa_h;
                        bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]) / kappa_h;
                        break;
                    case mode_tmy:
                        bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]) / kappa_h;
                        bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]) / kappa_h;
                        break;
                    case mode_tmz:
                        bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]) / kappa_h;
                        by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]) / kappa_h;
                        break;
                    case mode_tex:
                        bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]) / kappa_h;
                        bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]) / kappa_h;
                        break;
                    case mode_tey:
                        by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]) / kappa_h;
                        by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]) / kappa_h;
                        break;
                    case mode_tez:
                        bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]) / kappa_h;
                        bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]) / kappa_h;
                        break;
                }

                switch (mode)
                {
                    case mode_full:
                        bx[x][y][z] = bx[x][y][z] - d_t * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start];
                        bx[x][y][z] = bx[x][y][z] + d_t * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start];
                        by[x][y][z] = by[x][y][z] - d_t * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start];
                        by[x][y][z] = by[x][y][z] + d_t * psi[partition_index].exz[x-x_start][y-y_start][z-z_start];
                        bz[x][y][z] = bz[x][y][z] - d_t * psi[partition_index].exy[x-x_start][y-y_start][z-z_start];
                        bz[x][y][z] = bz[x][y][z] + d_t * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tmx:
                        by[x][y][z] = by[x][y][z] - d_t * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start];
                        bz[x][y][z] = bz[x][y][z] + d_t * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tmy:
                        bx[x][y][z] = bx[x][y][z] + d_t * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start];
                        bz[x][y][z] = bz[x][y][z] - d_t * psi[partition_index].exy[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tmz:
                        bx[x][y][z] = bx[x][y][z] - d_t * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start];
                        by[x][y][z] = by[x][y][z] + d_t * psi[partition_index].exz[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tex:
                        bx[x][y][z] = bx[x][y][z] - d_t * psi[partition_index].eyz[x-x_start][y-y_start][z-z_start];
                        bx[x][y][z] = bx[x][y][z] + d_t * psi[partition_index].ezy[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tey:
                        by[x][y][z] = by[x][y][z] - d_t * psi[partition_index].ezx[x-x_start][y-y_start][z-z_start];
                        by[x][y][z] = by[x][y][z] + d_t * psi[partition_index].exz[x-x_start][y-y_start][z-z_start];
                        break;
                    case mode_tez:
                        bz[x][y][z] = bz[x][y][z] - d_t * psi[partition_index].exy[x-x_start][y-y_start][z-z_start];
                        bz[x][y][z] = bz[x][y][z] + d_t * psi[partition_index].eyx[x-x_start][y-y_start][z-z_start];
                        break;
                }
            }
        }
    }
}

