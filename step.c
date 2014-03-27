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

static void update_b();
static void update_d();
static void b2h();
static void d2e();

void get_h()
{
    update_b();
    if (pml_type == pml_cpml) cpml_get_b();
    b2h();
    if (pml_type == pml_pml) pml_get_h();
}

void get_e()
{
    if (classical == 1)
        classical_e();
    else
    {
        update_d();
        if (pml_type == pml_cpml) 
            cpml_get_d();
        d2e();
    }
    if (pml_type == pml_pml) 
        pml_get_e();
}

static void update_b()
{
    for (z = 0; z < total_z; z++) 
        for (y = 0; y < total_y; y++) 
            for (x = 0; x < total_x; x++) 
            {
                if (!in_partition_main(x, y, z))
                    continue;

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

static void update_d()
{
    for (z = 0; z < total_z; z++) 
        for (y = 0; y < total_y; y++) 
            for (x = 0; x < total_x; x++) 
            {
                if (!in_partition_main(x, y, z))
                    continue;

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

void apply_pmc()
{
    if (partition_data[partition_z0].boundary_type == boundary_pmc)
        for (x = 0; x < total_x; x++)
            for (y = 0; y < total_y; y++)
            {
                if (ex)
                    ex[x][y][0] = (ex[x][y][1] * 4.0 - ex[x][y][2]) / 3.0;
                if (ey)
                    ey[x][y][0] = (ey[x][y][1] * 4.0 - ey[x][y][2]) / 3.0;
                if (ez)
                    ez[x][y][0] = (ez[x][y][1] * 4.0 - ez[x][y][2]) / 3.0;
            }
    if (partition_data[partition_z1].boundary_type == boundary_pmc)
        for (x = 0; x < total_x; x++)
            for (y = 0; y < total_y; y++)
            {
                if (ex)
                    ex[x][y][total_z] = (ex[x][y][total_z-1] * 4.0 - ex[x][y][total_z-2]) / 3.0;
                if (ey)
                    ey[x][y][total_z] = (ey[x][y][total_z-1] * 4.0 - ey[x][y][total_z-2]) / 3.0;
                if (ez)
                    ez[x][y][total_z] = (ez[x][y][total_z-1] * 4.0 - ez[x][y][total_z-2]) / 3.0;
            }
    if (partition_data[partition_x0].boundary_type == boundary_pmc)
        for (y = 0; y < total_y; y++)
            for (z = 0; z < total_z; z++)
            {
                if (ex)
                    ex[0][y][z] = (ex[1][y][z] * 4.0 - ex[2][y][z]) / 3.0;
                if (ey)
                    ey[0][y][z] = (ey[1][y][z] * 4.0 - ey[2][y][z]) / 3.0;
                if (ez)
                    ez[0][y][z] = (ez[1][y][z] * 4.0 - ez[2][y][z]) / 3.0;
            }
    if (partition_data[partition_x1].boundary_type == boundary_pmc)
        for (y = 0; y < total_y; y++)
            for (z = 0; z < total_z; z++)
            {
                if (ex)
                    ex[total_x][y][z] = (ex[total_x - 1][y][z] * 4.0 - ex[total_x - 2][y][z]) / 3.0;
                if (ey)
                    ey[total_x][y][z] = (ey[total_x - 1][y][z] * 4.0 - ey[total_x - 2][y][z]) / 3.0;
                if (ez)
                    ez[total_x][y][z] = (ez[total_x - 1][y][z] * 4.0 - ez[total_x - 2][y][z]) / 3.0;
            }
    if (partition_data[partition_y0].boundary_type == boundary_pmc && total_y > 2)
        for (x = 0; x < total_x; x++)
            for (z = 0; z < total_z; z++)
            {
                if (ex)
                    ex[x][0][z] = (ex[x][1][z] * 4.0 - ex[x][2][z]) / 3.0;
                if (ey)
                    ey[x][0][z] = (ey[x][1][z] * 4.0 - ey[x][2][z]) / 3.0;
                if (ez)
                    ez[x][0][z] = (ez[x][1][z] * 4.0 - ez[x][2][z]) / 3.0;
            }
    if (partition_data[partition_y1].boundary_type == boundary_pmc && total_y > 2)
        for (x = 0; x < total_x; x++)
            for (z = 0; z < total_z; z++)
            {
                if (ex)
                    ex[x][total_y][z] = (ex[x][total_y-1][z] * 4.0 - ex[x][total_y-2][z]) / 3.0;
                if (ey)
                    ey[x][total_y][z] = (ey[x][total_y-1][z] * 4.0 - ey[x][total_y-2][z]) / 3.0;
                if (ez)
                    ez[x][total_y][z] = (ez[x][total_y-1][z] * 4.0 - ez[x][total_y-2][z]) / 3.0;
            }
}

void b2h()
{
    for (z = 0; z < total_z; z++) 
        for (y = 0; y < total_y; y++) 
            for (x = 0; x < total_x; x++) 
            {
                if (hx)
                    bx[x][y][z] = bx[x][y][z] - d_t * dipole_hx[x][y][z];
                if (hy)
                    by[x][y][z] = by[x][y][z] - d_t * dipole_hy[x][y][z];
                if (hz)
                    bz[x][y][z] = bz[x][y][z] - d_t * dipole_hz[x][y][z];

                if(pml_type == pml_pml && !in_partition_main(x,y,z))
                    continue;

                if (hx)
                    hx[x][y][z] = bx[x][y][z] / MU0;
                if (hy)
                    hy[x][y][z] = by[x][y][z] / MU0;
                if (hz)
                    hz[x][y][z] = bz[x][y][z] / MU0;
            }
}

void d2e()
{
    for (x = 0; x < total_x; x++) 
        for (y = 0; y < total_y; y++) 
            for (z = 0; z < total_z; z++) 
            {
                if(pml_type == pml_pml && !in_partition_main(x,y,z))
                    continue;

                if (is_dispersive(x, y, z))
                {
                    if (ex && (y != 0 || z != 0)) 
                        ex[x][y][z] = get_ade_ex(x, y, z, dx[x][y][z], ex[x][y][z]);
                    if (ey && (x != 0 || z != 0))
                        ey[x][y][z] = get_ade_ey(x, y, z, dy[x][y][z], ey[x][y][z]);
                    if (ez && (x != 0 || y != 0))
                        ez[x][y][z] = get_ade_ez(x, y, z, dz[x][y][z], ez[x][y][z]);
                    continue;
                }

                if (ex && (y != 0 || z != 0)) 
                {
                    dx[x][y][z] = dx[x][y][z] - d_t * dipole_ex[x][y][z];
                    ex[x][y][z] = dx[x][y][z]  / epsilon[x][y][z];
                }
                if (ey && (x != 0 || z != 0))
                {
                    dy[x][y][z] = dy[x][y][z] - d_t * dipole_ey[x][y][z];
                    ey[x][y][z] = dy[x][y][z]  / epsilon[x][y][z];
                }
                if (ez && (x != 0 || y != 0))
                {
                    dz[x][y][z] = dz[x][y][z] - d_t * dipole_ez[x][y][z];
                    ez[x][y][z] = dz[x][y][z]  / epsilon[x][y][z];
                }
            }
}

