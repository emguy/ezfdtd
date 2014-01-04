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
 *     File Name : excitation.c
 * Last Modified : Fri 09 Nov 2012 04:20:04 PM EST
 */

#include "tools.h"
#include "h5io.h"
#include "excitation.h"

int **source_points;
int *source_polarizations;
double **source_signals;

int **hard_points;
int *hard_polarizations;
double **hard_signals;

int total_input_ports;
int total_hards;
int total_point_sources;

InputPort input_ports[100];

int starting_index;

int setup_input_ports(char* file_name)
{
    int port_index;
    char attr_name[45];
    int status;
    int point_index;

    status = h5_get_attr(file_name, "excitations", "number_of_ports", &total_input_ports);
    inspect(status, "fail to get h5 attributes");
    if (total_input_ports == 0) return 1;
    status = h5_get_attr(file_name, "settings", "starting_index", &starting_index);
    inspect(status, "fail to get h5 attributes");

    for (port_index = 1; port_index <= total_input_ports; port_index++)
    {
        sprintf(attr_name, "port_%02d_number_of_points", port_index);
        status = h5_get_attr(file_name, "excitations", attr_name, &(input_ports[port_index].total_points));
        inspect(status, "fail to get h5 attributes");
        sprintf(attr_name, "port_%02d_polarization", port_index);
        status = h5_get_attr(file_name, "excitations", attr_name, &(input_ports[port_index].polarization));
        inspect(status, "fail to get h5 attributes");
        inspect(mode == 0 || input_ports[port_index].polarization == mode, "polarization does not match with current simulation mode");

        sprintf(attr_name, "/excitations/port_%02d_points", port_index);
        input_ports[port_index].source_points = (int **)h5_load2(file_name, attr_name, input_ports[port_index].total_points, 3);
        inspect(input_ports[port_index].source_points, "fail to load hdf5 dataset");

        for (point_index = 0; point_index < input_ports[port_index].total_points; point_index ++)
        {
            input_ports[port_index].source_points[point_index][0] = input_ports[port_index].source_points[point_index][0] - starting_index;
            input_ports[port_index].source_points[point_index][1] = input_ports[port_index].source_points[point_index][1] - starting_index;
            input_ports[port_index].source_points[point_index][2] = input_ports[port_index].source_points[point_index][2] - starting_index;
        }

        sprintf(attr_name, "/excitations/port_%02d_signal", port_index);
        input_ports[port_index].source_signal = (double *)h5_load1(file_name, attr_name, total_timesteps);
        inspect(input_ports[port_index].source_signal, "fail to load hdf5 dataset");
        sprintf(attr_name, "/excitations/port_%02d_distribution", port_index);
        input_ports[port_index].source_distribution = (double *)h5_load1(file_name, attr_name, input_ports[port_index].total_points);
        inspect(input_ports[port_index].source_distribution, "fail to load hdf5 dataset");
    }
    return 1;
}

int setup_point_sources(char* file_name)
{
    int status;
    int point_index;

    status = h5_get_attr(file_name, "excitations", "number_of_point_sources", &total_point_sources);
    inspect(status, "fail to load hdf5 dataset");
    if (total_point_sources == 0) return 1;

    printf("-----------\n");
    status = h5_get_attr(file_name, "settings", "starting_index", &starting_index);
    inspect(status, "fail to get h5 attributes");

    source_points = (int **)h5_load2(file_name, "/excitations/source_points", total_point_sources, 3);
    inspect(source_points, "fail to load hdf5 dataset");

    for (point_index = 0; point_index < total_point_sources; point_index ++)
    {
        source_points[point_index][0] = source_points[point_index][0] - starting_index;
        source_points[point_index][1] = source_points[point_index][1] - starting_index;
        source_points[point_index][2] = source_points[point_index][2] - starting_index;
    }
    source_polarizations = (int *)h5_load1(file_name, "/excitations/source_polarizations", total_point_sources);
    inspect(source_polarizations, "fail to load hdf5 dataset");
    source_signals = (double **)h5_load2(file_name, "/excitations/source_signals", total_timesteps, total_point_sources);
    inspect(source_signals, "fail to load hdf5 dataset");

    return 1;
}

int setup_hards(char* file_name)
{
    int status;
    int point_index;

    status = h5_get_attr(file_name, "excitations", "number_of_hards", &total_hards);
    inspect(status, "fail to load hdf5 dataset");
    if (total_hards == 0) return 1;
    status = h5_get_attr(file_name, "settings", "starting_index", &starting_index);
    inspect(status, "fail to get h5 attributes");

    hard_points = (int **)h5_load2(file_name, "/excitations/hard_points", total_hards, 3);
    inspect(hard_points, "fail to load hdf5 dataset");
    for (point_index = 0; point_index < total_hards; point_index ++)
    {
        hard_points[point_index][0] = hard_points[point_index][0] - starting_index;
        hard_points[point_index][1] = hard_points[point_index][1] - starting_index;
        hard_points[point_index][2] = hard_points[point_index][2] - starting_index;
    }
    hard_polarizations = (int *)h5_load1(file_name, "/excitations/hard_polarizations", total_hards);
    inspect(hard_polarizations, "fail to load hdf5 dataset");
    hard_signals = (double **)h5_load2(file_name, "/excitations/hard_signals", total_timesteps, total_hards);
    inspect(hard_signals, "fail to load hdf5 dataset");
    return 1;
}

int excite(int time_index)
{
    int x_main, y_main, z_main;
    int source_index;
    int port_index;

    for (port_index = 1; port_index <= total_input_ports; port_index++)
    {
        for (source_index = 0; source_index < input_ports[port_index].total_points; source_index++)
        {
            x_main = input_ports[port_index].source_points[source_index][0];
            y_main = input_ports[port_index].source_points[source_index][1];
            z_main = input_ports[port_index].source_points[source_index][2];

            if (input_ports[port_index].polarization == p_ex)
                dipole_ex[x_main][y_main][z_main] = input_ports[port_index].source_signal[time_index] 
                    * input_ports[port_index].source_distribution[source_index];
            else if (input_ports[port_index].polarization == p_ey)
                dipole_ey[x_main][y_main][z_main] = input_ports[port_index].source_signal[time_index] 
                    * input_ports[port_index].source_distribution[source_index];
            else if (input_ports[port_index].polarization == p_ez)
                dipole_ez[x_main][y_main][z_main] = input_ports[port_index].source_signal[time_index] 
                    * input_ports[port_index].source_distribution[source_index];
            else if (input_ports[port_index].polarization == p_hx)
                dipole_hx[x_main][y_main][z_main] = input_ports[port_index].source_signal[time_index] 
                    * input_ports[port_index].source_distribution[source_index];
            else if (input_ports[port_index].polarization == p_hy)
                dipole_hy[x_main][y_main][z_main] = input_ports[port_index].source_signal[time_index] 
                    * input_ports[port_index].source_distribution[source_index];
            else if (input_ports[port_index].polarization == p_hz)
                dipole_hz[x_main][y_main][z_main] = input_ports[port_index].source_signal[time_index] 
                    * input_ports[port_index].source_distribution[source_index];
        }
    }
    for (source_index = 0; source_index < total_point_sources; source_index++)
    {
        x_main = source_points[source_index][0];
        y_main = source_points[source_index][1];
        z_main = source_points[source_index][2];
        if (source_polarizations[source_index] == p_ex)
            dipole_ex[x_main][y_main][z_main] = source_signals[time_index][source_index];
        else if (source_polarizations[source_index] == p_ey)
            dipole_ey[x_main][y_main][z_main] = source_signals[time_index][source_index];
        else if (source_polarizations[source_index] == p_ez)
            dipole_ez[x_main][y_main][z_main] = source_signals[time_index][source_index];
        else if (source_polarizations[source_index] == p_hx)
            dipole_hx[x_main][y_main][z_main] = source_signals[time_index][source_index];
        else if (source_polarizations[source_index] == p_hy)
        {
            dipole_hy[x_main][y_main][z_main] = source_signals[time_index][source_index];
            //printf("[%d %d %d] %e\n", x_main, y_main, z_main, dipole_hy[x_main][y_main][z_main]);
        }
        else if (source_polarizations[source_index] == p_hz)
            dipole_hz[x_main][y_main][z_main] = source_signals[time_index][source_index];
    }
    return 1;
}

int apply_ehards(int time_index)
{
    int x, y, z;
    int hard_index;

    for (hard_index = 0; hard_index < total_hards; hard_index++)
    {
        x = hard_points[hard_index][0] + partition_data[partition_main].x_start;
        y = hard_points[hard_index][1] + partition_data[partition_main].y_start;
        z = hard_points[hard_index][2] + partition_data[partition_main].z_start;
        if (hard_polarizations[hard_index] == p_ex)
            ex[x][y][z] = hard_signals[time_index][hard_index];
        else if (hard_polarizations[hard_index] == p_ey)
            ey[x][y][z] = hard_signals[time_index][hard_index];
        else if (hard_polarizations[hard_index] == p_ez)
            ez[x][y][z] = hard_signals[time_index][hard_index];
    }
    //for (hard_index = 0; hard_index < abc_size; hard_index++)
    //{
    //        ey[100][hard_index][0] = 0;
    //}
    return 1;
}


int apply_hhards(int time_index)
{
    int x, y, z;
    int hard_index;

    for (hard_index = 0; hard_index < total_hards; hard_index++)
    {
        x = hard_points[hard_index][0] + partition_data[partition_main].x_start;
        y = hard_points[hard_index][1] + partition_data[partition_main].y_start;
        z = hard_points[hard_index][2] + partition_data[partition_main].z_start;
        if (hard_polarizations[hard_index] == p_hx)
            hx[x][y][z] = hard_signals[time_index][hard_index];
        else if (hard_polarizations[hard_index] == p_hy)
            hy[x][y][z] = hard_signals[time_index][hard_index];
        else if (hard_polarizations[hard_index] == p_hz)
            hz[x][y][z] = hard_signals[time_index][hard_index];
    }
    //for (hard_index = 0; hard_index < abc_size; hard_index++)
    //{
    //        ey[100][hard_index][0] = 0;
    //}
    return 1;
}
