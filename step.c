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
 *     File Name : step.c
 * Last Modified : Fri 07 Dec 2012 10:30:15 PM EST
 */
#include "tools.h"
#include "domain.h"
#include "ade.h"
#include "step.h"
#include "excitation.h"
#include "cpml.h"
#include "pml.h"
#include "classical.h"

int x, y, z;
int x_main, y_main, z_main;

static void update_b();
static void update_d();
static void b2h();
static void d2e();

void get_h()
{
    update_b();
    if (pml_type != 0) cpml_get_b();
    else pml_get_h();
    b2h();
}

void get_e()
{
    if (classical != 1) update_d();
    else classical_e();
    if (pml_type != 0) cpml_get_d();
    else pml_get_e();
    d2e();
}

static void update_b()
{
    for (z = partition_data[partition_main].z_start; z < partition_data[partition_main].z_stop; z++) 
    {
        for (y = partition_data[partition_main].y_start; y < partition_data[partition_main].y_stop; y++) 
        {
            for (x = partition_data[partition_main].x_start; x < partition_data[partition_main].x_stop; x++) 
            {
                switch (mode)
                {
                    case mode_full:
                        bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]);
                        bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]);
                        by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]);
                        by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]);
                        bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]);
                        bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]);
                        break;
                    case mode_tmx:
                        by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]);
                        bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]);
                        break;
                    case mode_tmy:
                        bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]);
                        bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]);
                        break;
                    case mode_tmz:
                        bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]);
                        by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]);
                        break;
                    case mode_tex:
                        bx[x][y][z] = bx[x][y][z] + d_tz * (ey[x][y][z+1] - ey[x][y][z]);
                        bx[x][y][z] = bx[x][y][z] - d_ty * (ez[x][y+1][z] - ez[x][y][z]);
                        break;
                    case mode_tey:
                        by[x][y][z] = by[x][y][z] + d_tx * (ez[x+1][y][z] - ez[x][y][z]);
                        by[x][y][z] = by[x][y][z] - d_tz * (ex[x][y][z+1] - ex[x][y][z]);
                        break;
                    case mode_tez:
                        bz[x][y][z] = bz[x][y][z] + d_ty * (ex[x][y+1][z] - ex[x][y][z]);
                        bz[x][y][z] = bz[x][y][z] - d_tx * (ey[x+1][y][z] - ey[x][y][z]);
                        break;
                }
            }
        }
    }
}

static void update_d()
{
    for (z = partition_data[partition_main].z_start; z < partition_data[partition_main].z_stop; z++) 
    {
        for (y = partition_data[partition_main].y_start; y < partition_data[partition_main].y_stop; y++) 
        {
            for (x = partition_data[partition_main].x_start; x < partition_data[partition_main].x_stop; x++) 
            {
                switch (mode)
                {
                    case mode_full:
                        if (y != 0)
                            dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]);
                        if (z != 0)
                            dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]);
                        if (z != 0)
                            dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]);
                        if (x != 0)
                            dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]);
                        if (x != 0)
                            dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]);
                        if (y != 0)
                            dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case mode_tmx:
                        if (y != 0)
                            dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]);
                        if (z != 0)
                            dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]);
                        break;
                    case mode_tmy:
                        if (z != 0)
                            dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]);
                        if (x != 0)
                            dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]);
                        break;
                    case mode_tmz:
                        if (x != 0)
                            dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]);
                        if (y != 0)
                            dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case mode_tex:
                        if (z != 0)
                            dy[x][y][z] = dy[x][y][z] + d_tz * (hx[x][y][z] - hx[x][y][z-1]);
                        if (y != 0)
                            dz[x][y][z] = dz[x][y][z] - d_ty * (hx[x][y][z] - hx[x][y-1][z]);
                        break;
                    case mode_tey:
                        if (z != 0)
                            dx[x][y][z] = dx[x][y][z] - d_tz * (hy[x][y][z] - hy[x][y][z-1]);
                        if (x != 0)
                            dz[x][y][z] = dz[x][y][z] + d_tx * (hy[x][y][z] - hy[x-1][y][z]);
                        break;
                    case mode_tez:
                        if (y != 0)
                            dx[x][y][z] = dx[x][y][z] + d_ty * (hz[x][y][z] - hz[x][y-1][z]);
                        if (x != 0)
                            dy[x][y][z] = dy[x][y][z] - d_tx * (hz[x][y][z] - hz[x-1][y][z]);
                        break;
                }
            }
        }
    }
}

void apply_pmc()
{
    if (partition_data[partition_z0].boundary_type == boundary_pmc)
        for (x = 0; x < total_length_x; x++)
            for (y = 0; y < total_length_y; y++)
            {
                if (ex)
                    ex[x][y][0] = (ex[x][y][1] * 4.0 - ex[x][y][2]) / 3.0;
                if (ey)
                    ey[x][y][0] = (ey[x][y][1] * 4.0 - ey[x][y][2]) / 3.0;
                if (ez)
                    ez[x][y][0] = (ez[x][y][1] * 4.0 - ez[x][y][2]) / 3.0;
            }
    if (partition_data[partition_z1].boundary_type == boundary_pmc)
        for (x = 0; x < total_length_x; x++)
            for (y = 0; y < total_length_y; y++)
            {
                if (ex)
                    ex[x][y][total_length_z] = (ex[x][y][total_length_z-1] * 4.0 - ex[x][y][total_length_z-2]) / 3.0;
                if (ey)
                    ey[x][y][total_length_z] = (ey[x][y][total_length_z-1] * 4.0 - ey[x][y][total_length_z-2]) / 3.0;
                if (ez)
                    ez[x][y][total_length_z] = (ez[x][y][total_length_z-1] * 4.0 - ez[x][y][total_length_z-2]) / 3.0;
            }
    if (partition_data[partition_x0].boundary_type == boundary_pmc)
        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                if (ex)
                    ex[0][y][z] = (ex[1][y][z] * 4.0 - ex[2][y][z]) / 3.0;
                if (ey)
                    ey[0][y][z] = (ey[1][y][z] * 4.0 - ey[2][y][z]) / 3.0;
                if (ez)
                    ez[0][y][z] = (ez[1][y][z] * 4.0 - ez[2][y][z]) / 3.0;
            }
    if (partition_data[partition_x1].boundary_type == boundary_pmc)
        for (y = 0; y < total_length_y; y++)
            for (z = 0; z < total_length_z; z++)
            {
                if (ex)
                    ex[total_length_x][y][z] = (ex[total_length_x - 1][y][z] * 4.0 - ex[total_length_x - 2][y][z]) / 3.0;
                if (ey)
                    ey[total_length_x][y][z] = (ey[total_length_x - 1][y][z] * 4.0 - ey[total_length_x - 2][y][z]) / 3.0;
                if (ez)
                    ez[total_length_x][y][z] = (ez[total_length_x - 1][y][z] * 4.0 - ez[total_length_x - 2][y][z]) / 3.0;
            }
    if (partition_data[partition_y0].boundary_type == boundary_pmc && total_length_y > 2)
        for (x = 0; x < total_length_x; x++)
            for (z = 0; z < total_length_z; z++)
            {
                if (ex)
                    ex[x][0][z] = (ex[x][1][z] * 4.0 - ex[x][2][z]) / 3.0;
                if (ey)
                    ey[x][0][z] = (ey[x][1][z] * 4.0 - ey[x][2][z]) / 3.0;
                if (ez)
                    ez[x][0][z] = (ez[x][1][z] * 4.0 - ez[x][2][z]) / 3.0;
            }
    if (partition_data[partition_y1].boundary_type == boundary_pmc && total_length_y > 2)
        for (x = 0; x < total_length_x; x++)
            for (z = 0; z < total_length_z; z++)
            {
                if (ex)
                    ex[x][total_length_y][z] = (ex[x][total_length_y-1][z] * 4.0 - ex[x][total_length_y-2][z]) / 3.0;
                if (ey)
                    ey[x][total_length_y][z] = (ey[x][total_length_y-1][z] * 4.0 - ey[x][total_length_y-2][z]) / 3.0;
                if (ez)
                    ez[x][total_length_y][z] = (ez[x][total_length_y-1][z] * 4.0 - ez[x][total_length_y-2][z]) / 3.0;
            }
}

void b2h()
{
    for (z = 0; z < total_length_z; z++) 
    {
        for (y = 0; y < total_length_y; y++) 
        {
            for (x = 0; x < total_length_x; x++) 
            {
                if (get_partition(x, y, z) == partition_main)
                {
                    x_main = x - partition_data[partition_main].x_start;
                    y_main = y - partition_data[partition_main].y_start;
                    z_main = z - partition_data[partition_main].z_start;

                    if (hx)
                        bx[x][y][z] = bx[x][y][z] - d_t * dipole_hx[x_main][y_main][z_main];
                    if (hy)
                        by[x][y][z] = by[x][y][z] - d_t * dipole_hy[x_main][y_main][z_main];
                    if (hz)
                        bz[x][y][z] = bz[x][y][z] - d_t * dipole_hz[x_main][y_main][z_main];
                    if (hx)
                        hx[x][y][z] = bx[x][y][z] / MU0;
                    if (hy)
                    {
                        hy[x][y][z] = by[x][y][z] / MU0;
                    }
                    if (hz)
                        hz[x][y][z] = bz[x][y][z] / MU0;
                }
                else if (pml_type != 0)
                {
                    if (hx)
                        hx[x][y][z] = bx[x][y][z] / MU0;
                    if (hy)
                        hy[x][y][z] = by[x][y][z] / MU0;
                    if (hz)
                        hz[x][y][z] = bz[x][y][z] / MU0;
                }
            }
        }
    }
}

void d2e()
{
    for (x = 0; x < total_length_x; x++) 
    {
        for (y = 0; y < total_length_y; y++) 
        {
            for (z = 0; z < total_length_z; z++) 
            {
                if (classical != 1 && get_partition(x, y, z) == 0)
                {
                    x_main = x - partition_data[partition_main].x_start;
                    y_main = y - partition_data[partition_main].y_start;
                    z_main = z - partition_data[partition_main].z_start;

                    if (!is_dispersive(x_main, y_main, z_main))
                    {
                        if (ex)
                        {
                            dx[x][y][z] = dx[x][y][z] - d_t * dipole_ex[x_main][y_main][z_main];
                            ex[x][y][z] = dx[x][y][z]  / epsilon[x_main][y_main][z_main];
                        }
                        if (ey)
                        {
                            dy[x][y][z] = dy[x][y][z] - d_t * dipole_ey[x_main][y_main][z_main];
                            ey[x][y][z] = dy[x][y][z]  / epsilon[x_main][y_main][z_main];
                        }
                        if (ez)
                        {
                            dz[x][y][z] = dz[x][y][z] - d_t * dipole_ez[x_main][y_main][z_main];
                            ez[x][y][z] = dz[x][y][z]  / epsilon[x_main][y_main][z_main];
                        }
                    }
                    else
                    {
                        if (ex)
                            ex[x][y][z] = get_ade_ex(x_main, y_main, z_main, dx[x][y][z], ex[x][y][z]);
                        if (ey)
                            ey[x][y][z] = get_ade_ey(x_main, y_main, z_main, dy[x][y][z], ey[x][y][z]);
                        if (ez)
                            ez[x][y][z] = get_ade_ez(x_main, y_main, z_main, dz[x][y][z], ez[x][y][z]);
                    }
                }
                else if (pml_type != 0 && get_partition(x, y, z) != 0)
                {
                    if (ex)
                        ex[x][y][z] = dx[x][y][z] / EPSILON0;
                    if (ey)
                        ey[x][y][z] = dy[x][y][z] / EPSILON0;
                    if (ez)
                        ez[x][y][z] = dz[x][y][z] / EPSILON0;
                }
            }
        }
    }
}
