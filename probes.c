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
 * Last Modified : Fri 21 Mar 2014 12:17:34 AM EDT
 */

#include "tools.h"
#include "ezfdtd.h"
#include "h5io.h"
#include "mem.h"
#include "domain.h"
#include "probes.h"

OutputPort output_ports[OUTPUTPORTS_MAX];
FieldPlane field_planes[10];

int x, y, z;
int starting_index;

int total_output_ports;
int total_planes;

int setup_planes(char* file_name, char* output_file_name)
{
    char attr_name[30];
    int plane_index;
    int status;
    int length_x, length_y, length_z;

    status = h5_get_attr(file_name, "outputs", "number_of_field_planes", &total_planes);
    inspect(status, "fail to get h5 attributes");
    if (total_planes == 0) return 1;
    status = h5_get_attr(file_name, "settings", "starting_index", &starting_index);
    inspect(status, "fail to get h5 attributes");

    for (plane_index = 1; plane_index <= total_planes; plane_index++)
    {
        sprintf(attr_name, "plane_%02d_polarization", plane_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].polarization));
        inspect(status, "fail to get h5 attributes");

        sprintf(attr_name, "plane_%02d_dim", plane_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].dim));
        inspect(status, "fail to get h5 attributes");

        sprintf(attr_name, "plane_%02d_slice_index", plane_index);
        status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].slice_index));
        inspect(status, "fail to get h5 attributes");

        field_planes[plane_index].slice_index -= starting_index;

        if (field_planes[plane_index].dim == dim_x)
        {
            sprintf(attr_name, "plane_%02d_y_start", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].y_start));
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].y_start -= starting_index;

            sprintf(attr_name, "plane_%02d_z_start", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].z_start));
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].z_start -= starting_index;

            sprintf(attr_name, "plane_%02d_length_y", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &length_y);
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].length_y = length_y;
            field_planes[plane_index].y_stop = field_planes[plane_index].y_start + length_y;

            sprintf(attr_name, "plane_%02d_length_z", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &length_z);
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].length_z = length_z;
            field_planes[plane_index].z_stop = field_planes[plane_index].z_start + length_z;

            sprintf(attr_name, "plane_%02d", plane_index);
            field_planes[plane_index].field_buffer = (double **)contiguous_mem2(type_double, length_y, length_z);
            status = h5_empty3(output_file_name, attr_name, type_double, length_y, length_z, total_timesteps);
            inspect(status, "fail to create an empty dset");
        }
        else if (field_planes[plane_index].dim == dim_y)
        {
            sprintf(attr_name, "plane_%02d_z_start", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].z_start));
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].z_start -= starting_index;

            sprintf(attr_name, "plane_%02d_x_start", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].x_start));
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].x_start -= starting_index;

            sprintf(attr_name, "plane_%02d_length_z", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &length_z);
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].length_z = length_z;
            field_planes[plane_index].z_stop = field_planes[plane_index].z_start + length_z;

            sprintf(attr_name, "plane_%02d_length_x", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &length_x);
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].length_x = length_x;
            field_planes[plane_index].x_stop = field_planes[plane_index].x_start + length_x;

            sprintf(attr_name, "plane_%02d", plane_index);
            field_planes[plane_index].field_buffer = (double **)contiguous_mem2(type_double, length_x, length_z);
            status = h5_empty3(output_file_name, attr_name, type_double, length_x, length_z, total_timesteps);
            inspect(status, "fail to create an empty dset");
        }
        else if (field_planes[plane_index].dim == dim_z)
        {
            sprintf(attr_name, "plane_%02d_x_start", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].x_start));
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].x_start -= starting_index;

            sprintf(attr_name, "plane_%02d_y_start", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &(field_planes[plane_index].y_start));
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].y_start -= starting_index;

            sprintf(attr_name, "plane_%02d_length_x", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &length_x);
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].length_x = length_x;
            field_planes[plane_index].x_stop = field_planes[plane_index].x_start + length_x;

            sprintf(attr_name, "plane_%02d_length_y", plane_index);
            status = h5_get_attr(file_name, "outputs", attr_name, &length_y);
            inspect(status, "fail to get h5 attributes");
            field_planes[plane_index].length_y = length_y;
            field_planes[plane_index].y_stop = field_planes[plane_index].y_start + length_y;

            sprintf(attr_name, "plane_%02d", plane_index);
            field_planes[plane_index].field_buffer = (double **)contiguous_mem2(type_double, length_x, length_y);
            status = h5_empty3(output_file_name, attr_name, type_double, length_x, length_y, total_timesteps);
            inspect(status, "fail to create an empty dset");
        }
    }
    return 1;
}

int update_planes(char* file_name, int time_index)
{
    int plane_index;
    char attr_name[30];
    int status;
    int x_start, y_start, z_start;
    int x_stop, y_stop, z_stop;
    int length_x, length_y, length_z;

    if (total_planes == 0) return 1;
    for (plane_index = 1; plane_index <= total_planes; plane_index++)
    {
        if (field_planes[plane_index].dim == 1)
        {
            y_start = field_planes[plane_index].y_start;
            z_start = field_planes[plane_index].z_start;
            y_stop = field_planes[plane_index].y_stop;
            z_stop = field_planes[plane_index].z_stop;

            x = field_planes[plane_index].slice_index;

            for (z = z_start; z < z_stop; z++)
                for (y = y_start; y < y_stop; y++)
                {
                    switch (field_planes[plane_index].polarization)
                    {
                        case (p_ex):
                            field_planes[plane_index].field_buffer[y - y_start][z - z_start] = ex[x][y][z];
                            break;
                        case (p_ey):
                            field_planes[plane_index].field_buffer[y - y_start][z - z_start] = ey[x][y][z];
                            break;
                        case (p_ez):
                            field_planes[plane_index].field_buffer[y - y_start][z - z_start] = ez[x][y][z];
                            break;
                        case (p_hx):
                            field_planes[plane_index].field_buffer[y - y_start][z - z_start] = hx[x][y][z];
                            break;
                        case (p_hy):
                            field_planes[plane_index].field_buffer[y - y_start][z - z_start] = hy[x][y][z];
                            break;
                        case (p_hz):
                            field_planes[plane_index].field_buffer[y - y_start][z - z_start] = hz[x][y][z];
                            break;
                    }
                }
        }

        if (field_planes[plane_index].dim == 2)
        {
            x_start = field_planes[plane_index].x_start;
            z_start = field_planes[plane_index].z_start;
            x_stop = field_planes[plane_index].x_stop;
            z_stop = field_planes[plane_index].z_stop;

            y = field_planes[plane_index].slice_index;

            for (z = z_start; z < z_stop; z++)
                for (x = x_start; x < x_stop; x++)
                {
                    switch (field_planes[plane_index].polarization)
                    {
                        case (p_ex):
                            field_planes[plane_index].field_buffer[x - x_start][z - z_start] = ex[x][y][z];
                            break;
                        case (p_ey):
                            field_planes[plane_index].field_buffer[x - x_start][z - z_start] = ey[x][y][z];
                            break;
                        case (p_ez):
                            field_planes[plane_index].field_buffer[x - x_start][z - z_start] = ez[x][y][z];
                            break;
                        case (p_hx):
                            field_planes[plane_index].field_buffer[x - x_start][z - z_start] = hx[x][y][z];
                            break;
                        case (p_hy):
                            field_planes[plane_index].field_buffer[x - x_start][z - z_start] = hy[x][y][z];
                            break;
                        case (p_hz):
                            field_planes[plane_index].field_buffer[x - x_start][z - z_start] = hz[x][y][z];
                            break;
                    }
                }
        }

        if (field_planes[plane_index].dim == 3)
        {
            x_start = field_planes[plane_index].x_start;
            y_start = field_planes[plane_index].y_start;
            x_stop = field_planes[plane_index].x_stop;
            y_stop = field_planes[plane_index].y_stop;

            z = field_planes[plane_index].slice_index;

            for (x = x_start; x < x_stop; x++)
                for (y = y_start; y < y_stop; y++)
                {
                    switch (field_planes[plane_index].polarization)
                    {
                        case (p_ex):
                            field_planes[plane_index].field_buffer[x - x_start][y - y_start] = ex[x][y][z];
                            break;
                        case (p_ey):
                            field_planes[plane_index].field_buffer[x - x_start][y - y_start] = ey[x][y][z];
                            break;
                        case (p_ez):
                            field_planes[plane_index].field_buffer[x - x_start][y - y_start] = ez[x][y][z];
                            break;
                        case (p_hx):
                            field_planes[plane_index].field_buffer[x - x_start][y - y_start] = hx[x][y][z];
                            break;
                        case (p_hy):
                            field_planes[plane_index].field_buffer[x - x_start][y - y_start] = hy[x][y][z];
                            break;
                        case (p_hz):
                            field_planes[plane_index].field_buffer[x - x_start][y - y_start] = hz[x][y][z];
                            break;
                    }
                }
        }

        for (plane_index = 1; plane_index <= total_planes; plane_index++)
        {
            length_x = field_planes[plane_index].length_x;
            length_y = field_planes[plane_index].length_y;
            length_z = field_planes[plane_index].length_z;

            sprintf(attr_name, "plane_%02d", plane_index);
            switch (field_planes[plane_index].dim)
            {
                case (dim_x):
                    status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, length_y, length_z, time_index);
                    inspect(status, "fail to write h5 hyperslab");
                    break;
                case (dim_y):
                    status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, length_x, length_z, time_index);
                    inspect(status, "fail to write h5 hyperslab");
                    break;
                case (dim_z):
                    status = h5_slab(file_name, attr_name, type_double, (void **)field_planes[plane_index].field_buffer, length_x, length_y, time_index);
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

    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        for (probe_index = 0; probe_index < output_ports[port_index].total_points; probe_index++)
        {
            x = output_ports[port_index].probe_points[probe_index][0];
            y = output_ports[port_index].probe_points[probe_index][1];
            z = output_ports[port_index].probe_points[probe_index][2];
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
