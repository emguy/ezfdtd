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
 *     File Name : dft.c
 * Last Modified : Thu 11 Oct 2012 04:53:30 PM EDT
 */

#include "tools.h"
#include "math.h"
#include "mem.h"
#include "h5io.h"
#include "dft.h"
#include "probes.h"

int total_frequencies;
double  *frequencies;

double flag_dft;

int setup_dft(char* file_name)
{
    int port_index;
    int total_points;
    int status;

    status = h5_get_attr(file_name, "DFT", "enable_DFT", &flag_dft);
    inspect(status, "fail to get h5 attributes");

    if (flag_dft == 0) return 1;
    status = h5_get_attr(file_name, "DFT", "number_of_frequencies", &total_frequencies);
    inspect(status, "fail to get h5 attributes");


    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        total_points = output_ports[port_index].total_points;
        frequencies = (double *)h5_load1(file_name, "/DFT/frequencies", total_frequencies);
        inspect(frequencies, "fail to load hdf5 dataset");
        output_ports[port_index].dft_real = (double **)mem2(type_double, total_frequencies, total_points);
        inspect(output_ports[port_index].dft_real, "fail to allocate memory");
        output_ports[port_index].dft_imag = (double **)mem2(type_double, total_frequencies, total_points);
        inspect(output_ports[port_index].dft_imag, "fail to allocate memory");
    }
    return 1;
}

int update_dft(int time_index)
{
    int port_index;
    int probe_index;
    int total_points;
    int frequency_index;
    double value;

    if (flag_dft == 0) return 1;
    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        total_points = output_ports[port_index].total_points;
        for (probe_index = 0; probe_index < total_points; probe_index++)
        {
            value = output_ports[port_index].fields[time_index][probe_index];
            for (frequency_index = 0; frequency_index < total_frequencies; frequency_index ++)
            {
                output_ports[port_index].dft_real[frequency_index][probe_index] 
                    = output_ports[port_index].dft_real[frequency_index][probe_index] 
                    + value *  cos(2.0 * PI * frequencies[frequency_index] * time_index * d_t);
                output_ports[port_index].dft_imag[frequency_index][probe_index] 
                    = output_ports[port_index].dft_imag[frequency_index][probe_index] 
                    - value *  sin(2.0 * PI * frequencies[frequency_index] * time_index * d_t);
            }
        }
    }
    return 1;
}

int write_dft(char* file_name)
{
    int port_index;
    char attr_name[30];
    int status;

    if (flag_dft == 0) return 1;
    for (port_index = 1; port_index <= total_output_ports; port_index++)
    {
        sprintf(attr_name, "port_%02d_spectrum_real", port_index);
        status = h5_write2(file_name, attr_name, type_double, (void **)output_ports[port_index].dft_real, total_frequencies, output_ports[port_index].total_points);
        inspect(status, "fail to write h5 dset");
        sprintf(attr_name, "port_%02d_spectrum_imag", port_index);
        status = h5_write2(file_name, attr_name, type_double, (void **)output_ports[port_index].dft_imag, total_frequencies, output_ports[port_index].total_points);
        inspect(status, "fail to write h5 dset");
    }
    return 1;
}

