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
 *     File Name : probes.c
 * Last Modified : Fri 12 Oct 2012 01:08:28 PM EDT
 */

#include "tools.h"
#include <string.h>
#include "ezfdtd.h"
#include "h5io.h"
#include "mem.h"
#include "domain.h"
#include "probes.h"

OutputPort output_ports[100];
FieldPlane field_planes[10];

int starting_index;

int total_output_ports;
int total_planes;

int setup_planes(char* file_name, char* output_file_name, int with_pml)
{
    char attr_name[30];
    int plane_index;
    int status;

    int total_x;
    int total_y;
    int total_z;

    total_x = total_length_x;
    total_y = total_length_y;
    total_z = total_length_z;
    if(!with_pml)
    {
        total_x = main_length_x;
        total_y = main_length_y;
        total_z = main_length_z;
    }

    status = h5_get_attr(file_name, "outputs", "number_of_field_planes", &total_planes);
    inspect(status, "fail to get h5 attributes");
    if (total_planes == 0) return 1;
    status = h5_get_attr(file_name, "settings", "starting_index", &starting_index);
    inspect(status, "fail to get h5 attributes");

    for (plane_index = 1; plane_index <= total_planes; plane_index++)
    {
        field_planes[plane_index].with_pml = with_pml;
        sprintf(attr_name, "plane_%02d_polarization", plane_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].polarization));
        inspect(status, "fail to get h5 attributes");
        sprintf(attr_name, "plane_%02d_dim", plane_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].dim));
        inspect(status, "fail to get h5 attributes");
        sprintf(attr_name, "plane_%02d_slice_index", plane_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].slice_index));
        inspect(status, "fail to get h5 attributes");
        field_planes[plane_index].slice_index = field_planes[plane_index].slice_index - starting_index;

        sprintf(attr_name, "plane_%02d", plane_index);
        if (field_planes[plane_index].dim == 1)
        {
            field_planes[plane_index].field_buffer = (double **)contiguous_mem2(type_double, total_y, total_z);
            status = h5_empty3(output_file_name, attr_name, type_double, total_y, total_z, total_timesteps);
            inspect(status, "fail to create an empty dset");
        }
        else if (field_planes[plane_index].dim == 2)
        {
            field_planes[plane_index].field_buffer = (double **)contiguous_mem2(type_double, total_x, total_z);
            status = h5_empty3(output_file_name, attr_name, type_double, total_x, total_z, total_timesteps);
            inspect(status, "fail to create an empty dset");
        }
        else if (field_planes[plane_index].dim == 3)
        {
            field_planes[plane_index].field_buffer = (double **)contiguous_mem2(type_double, total_x, total_y);
            status = h5_empty3(output_file_name, attr_name, type_double, total_x, total_y, total_timesteps);
            inspect(status, "fail to create an empty dset");
        }
    }
    return 1;
}

int update_planes(char* file_name, int time_index)
{
    int plane_index;
    int x, y, z;
    int x_main, y_main, z_main;
    char attr_name[30];
    int status;

    if (total_planes == 0) return 1;
    for (plane_index = 1; plane_index <= total_planes; plane_index++)
    {
        for (x = 0; x < total_length_x; x++)
        {
            for (y = 0; y < total_length_y; y++)
            {
                for (z = 0; z < total_length_z; z++)
                {
                    x_main = x - partition_data[partition_main].x_start;
                    y_main = y - partition_data[partition_main].y_start;
                    z_main = z - partition_data[partition_main].z_start;

                    if(field_planes[plane_index].with_pml == 0 && get_partition(x, y, z) == 0)
                    {

                        if (field_planes[plane_index].dim == 1 && x_main == field_planes[plane_index].slice_index)
                        {
                            if (field_planes[plane_index].polarization == 1)
                                field_planes[plane_index].field_buffer[y_main][z_main] = ex[x][y][z];
                            else if (field_planes[plane_index].polarization == 2)
                                field_planes[plane_index].field_buffer[y_main][z_main] = ey[x][y][z];
                            else if (field_planes[plane_index].polarization == 3)
                                field_planes[plane_index].field_buffer[y_main][z_main] = ez[x][y][z];
                            else if (field_planes[plane_index].polarization == 4)
                                field_planes[plane_index].field_buffer[y_main][z_main] = hx[x][y][z];
                            else if (field_planes[plane_index].polarization == 5)
                                field_planes[plane_index].field_buffer[y_main][z_main] = hy[x][y][z];
                            else if (field_planes[plane_index].polarization == 6)
                                field_planes[plane_index].field_buffer[y_main][z_main] = hz[x][y][z];
                        }
                        else if (field_planes[plane_index].dim == 2 && y_main == field_planes[plane_index].slice_index)
                        {
                            if (field_planes[plane_index].polarization == 1)
                                field_planes[plane_index].field_buffer[x_main][z_main] = ex[x][y][z];
                            else if (field_planes[plane_index].polarization == 2)
                                field_planes[plane_index].field_buffer[x_main][z_main] = ey[x][y][z];
                            else if (field_planes[plane_index].polarization == 3)
                                field_planes[plane_index].field_buffer[x_main][z_main] = ez[x][y][z];
                            else if (field_planes[plane_index].polarization == 4)
                                field_planes[plane_index].field_buffer[x_main][z_main] = hx[x][y][z];
                            else if (field_planes[plane_index].polarization == 5)
                                field_planes[plane_index].field_buffer[x_main][z_main] = hy[x][y][z];
                            else if (field_planes[plane_index].polarization == 6)
                                field_planes[plane_index].field_buffer[x_main][z_main] = hz[x][y][z];
                        }
                        else if (field_planes[plane_index].dim == 3 && z_main == field_planes[plane_index].slice_index)
                        {
                            if (field_planes[plane_index].polarization == 1)
                                field_planes[plane_index].field_buffer[x_main][y_main] = ex[x][y][z];
                            else if (field_planes[plane_index].polarization == 2)
                                field_planes[plane_index].field_buffer[x_main][y_main] = ey[x][y][z];
                            else if (field_planes[plane_index].polarization == 3)
                                field_planes[plane_index].field_buffer[x_main][y_main] = ez[x][y][z];
                            else if (field_planes[plane_index].polarization == 4)
                                field_planes[plane_index].field_buffer[x_main][y_main] = hx[x][y][z];
                            else if (field_planes[plane_index].polarization == 5)
                                field_planes[plane_index].field_buffer[x_main][y_main] = hy[x][y][z];
                            else if (field_planes[plane_index].polarization == 6)
                                field_planes[plane_index].field_buffer[x_main][y_main] = hz[x][y][z];
                        }
                    }
                    if (field_planes[plane_index].with_pml == 1)
                    {
                        if (field_planes[plane_index].dim == 1 && x_main == field_planes[plane_index].slice_index)
                        {
                            if (field_planes[plane_index].polarization == 1)
                                field_planes[plane_index].field_buffer[y][z] = ex[x][y][z];
                            else if (field_planes[plane_index].polarization == 2)
                                field_planes[plane_index].field_buffer[y][z] = ey[x][y][z];
                            else if (field_planes[plane_index].polarization == 3)
                                field_planes[plane_index].field_buffer[y][z] = ez[x][y][z];
                            else if (field_planes[plane_index].polarization == 4)
                                field_planes[plane_index].field_buffer[y][z] = hx[x][y][z];
                            else if (field_planes[plane_index].polarization == 5)
                                field_planes[plane_index].field_buffer[y][z] = hy[x][y][z];
                            else if (field_planes[plane_index].polarization == 6)
                                field_planes[plane_index].field_buffer[y][z] = hz[x][y][z];
                        }
                        else if (field_planes[plane_index].dim == 2 && y_main == field_planes[plane_index].slice_index)
                        {
                            if (field_planes[plane_index].polarization == 1)
                                field_planes[plane_index].field_buffer[x][z] = ex[x][y][z];
                            else if (field_planes[plane_index].polarization == 2)
                                field_planes[plane_index].field_buffer[x][z] = ey[x][y][z];
                            else if (field_planes[plane_index].polarization == 3)
                                field_planes[plane_index].field_buffer[x][z] = ez[x][y][z];
                            else if (field_planes[plane_index].polarization == 4)
                                field_planes[plane_index].field_buffer[x][z] = hx[x][y][z];
                            else if (field_planes[plane_index].polarization == 5)
                                field_planes[plane_index].field_buffer[x][z] = hy[x][y][z];
                            else if (field_planes[plane_index].polarization == 6)
                                field_planes[plane_index].field_buffer[x][z] = hz[x][y][z];
                        }
                        else if (field_planes[plane_index].dim == 3 && z_main == field_planes[plane_index].slice_index)
                        {
                            if (field_planes[plane_index].polarization == 1)
                                field_planes[plane_index].field_buffer[x][y] = ex[x][y][z];
                            else if (field_planes[plane_index].polarization == 2)
                                field_planes[plane_index].field_buffer[x][y] = ey[x][y][z];
                            else if (field_planes[plane_index].polarization == 3)
                                field_planes[plane_index].field_buffer[x][y] = ez[x][y][z];
                            else if (field_planes[plane_index].polarization == 4)
                                field_planes[plane_index].field_buffer[x][y] = hx[x][y][z];
                            else if (field_planes[plane_index].polarization == 5)
                                field_planes[plane_index].field_buffer[x][y] = hy[x][y][z];
                            else if (field_planes[plane_index].polarization == 6)
                                field_planes[plane_index].field_buffer[x][y] = hz[x][y][z];
                        }
                    }
                }
            }
        }

        for (plane_index = 1; plane_index <= total_planes; plane_index++)
        {
            sprintf(attr_name, "plane_%02d", plane_index);
            switch (field_planes[plane_index].dim)
            {
                case (dim_x):
                    if (field_planes[plane_index].with_pml)
                        status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, total_length_y, total_length_z, time_index);
                    else
                        status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, main_length_y, main_length_z, time_index);
                    inspect(status, "fail to write h5 hyperslab");
                    break;
                case (dim_y):
                    if (field_planes[plane_index].with_pml)
                        status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, total_length_x, total_length_z, time_index);
                    else
                        status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, main_length_x, main_length_z, time_index);
                    inspect(status, "fail to write h5 hyperslab");
                    break;
                case (dim_z):
                    if (field_planes[plane_index].with_pml)
                        status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, total_length_x, total_length_y, time_index);
                    else
                        status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, main_length_x, main_length_y, time_index);
                    inspect(status, "fail to write h5 hyperslab");
                    break;
            }
        }
    }
    return 1;
}

int setup_output_ports(char* file_name)
{
    int port_index;
    char attr_name[30];
    int status;
    int point_index;

    status = h5_get_attr(file_name, "outputs", "number_of_ports", &total_output_ports);

    if (total_output_ports == 0) return 1;
    status = h5_get_attr(file_name, "settings", "starting_index", &starting_index);
    inspect(status, "fail to get h5 attributes");

    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        sprintf(attr_name, "port_%02d_write_timedomain", port_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(output_ports[port_index].write_timedomain));
        inspect(status, "fail to get h5 attributes");
        sprintf(attr_name, "port_%02d_number_of_points", port_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(output_ports[port_index].total_points));
        inspect(status, "fail to get h5 attributes");
        sprintf(attr_name, "port_%02d_polarization", port_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(output_ports[port_index].polarization));
        inspect(status, "fail to get h5 attributes");
        output_ports[port_index].fields = (double **)mem2(type_double, total_timesteps, output_ports[port_index].total_points);
        inspect(output_ports[port_index].fields, "fail to allocate memory for output probes");
        sprintf(attr_name, "/outputs/port_%02d_points", port_index);
        output_ports[port_index].probe_points = (int **)h5_load2(file_name, attr_name, output_ports[port_index].total_points, 3);
        inspect(output_ports[port_index].probe_points, "fail to allocate memory for port points");
        for (point_index = 0; point_index < output_ports[port_index].total_points; point_index ++)
        {
            output_ports[port_index].probe_points[point_index][0] = output_ports[port_index].probe_points[point_index][0] - starting_index;
            output_ports[port_index].probe_points[point_index][1] = output_ports[port_index].probe_points[point_index][1] - starting_index;
            output_ports[port_index].probe_points[point_index][2] = output_ports[port_index].probe_points[point_index][2] - starting_index;
        }
    }
    return 1;
}

int update_ports(int time_index)
{
    int port_index;
    int probe_index;
    int x, y, z;

    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        for (probe_index = 0; probe_index < output_ports[port_index].total_points; probe_index++)
        {
            x = output_ports[port_index].probe_points[probe_index][0] + partition_data[partition_main].x_start;
            y = output_ports[port_index].probe_points[probe_index][1] + partition_data[partition_main].y_start;
            z = output_ports[port_index].probe_points[probe_index][2] + partition_data[partition_main].z_start;
            if (output_ports[port_index].polarization == 1)
                switch (mode)
                {
                    case 5:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (ex[x][y][z] + ex[x][y][z+1]);
                        break;
                    case 6:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (ex[x][y][z] + ex[x][y+1][z]);
                        break;
                    default:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.25 * (ex[x][y][z] + ex[x][y+1][z] + ex[x][y][z+1] 
                                    + ex[x][y+1][z+1]);
                }
            else if (output_ports[port_index].polarization == 2)
                switch (mode)
                {
                    case 4:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (ey[x][y][z] + ey[x][y][z+1]);
                        break;
                    case 6:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (ey[x][y][z] + ey[x+1][y][z]);
                        break;
                    default:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.25 * (ey[x][y][z] + ey[x+1][y][z] + ey[x][y][z+1] 
                                    + ey[x+1][y][z+1]);
                }
            else if (output_ports[port_index].polarization == 3)
                switch (mode)
                {
                    case 4:
                        output_ports[port_index].fields[time_index][probe_index]
                            = 0.5 * (ez[x][y][z] + ez[x][y+1][z]);
                        break;
                    case 5:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (ez[x][y][z] + ez[x+1][y][z]);
                        break;
                    default:
                        output_ports[port_index].fields[time_index][probe_index] = ez[x][y][z];
                        //output_ports[port_index].fields[time_index][probe_index]
                        //    = 0.25 * (ez[x][y][z] + ez[x+1][y][z] + ez[x][y+1][z] 
                        //           + ez[x+1][y+1][z]);
                }
            else if (output_ports[port_index].polarization == 4)
                switch (mode)
                {
                    case 4:
                        output_ports[port_index].fields[time_index][probe_index] = hx[x][y][z];
                        break;
                    default:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (hx[x][y][z] + hx[x+1][y][z]);
                }
            else if (output_ports[port_index].polarization == 5)
                switch (mode)
                {
                    case 5:
                        output_ports[port_index].fields[time_index][probe_index] = hy[x][y][z];
                        break;
                    default:
                        output_ports[port_index].fields[time_index][probe_index] 
                        = 0.5 * (hy[x][y][z] + hy[x][y+1][z]);
                }
            else if (output_ports[port_index].polarization == 6)
                switch (mode)
                {
                    case 6:
                        output_ports[port_index].fields[time_index][probe_index] = hz[x][y][z];
                        break;
                    default:
                        output_ports[port_index].fields[time_index][probe_index] 
                            = 0.5 * (hz[x][y][z] + hz[x][y][z+1]);
                }
        }
    }
    return 1;
}

int write_ports(char* file_name)
{
    int port_index;
    char attr_name[30];
    int status;

    if (total_output_ports == 0) return 1;
    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        if (! output_ports[port_index].write_timedomain) continue;
        sprintf(attr_name, "port_%02d", port_index);
        status = h5_write2(file_name, attr_name, type_double, (void **)output_ports[port_index].fields, total_timesteps, output_ports[port_index].total_points);
        inspect(status, "fail to write h5 dset");
    }
    return 1;
}
