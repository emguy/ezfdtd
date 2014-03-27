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
 *     File Name : ezfdtd.c
 * Last Modified : Sun 07 Oct 2012 12:41:45 AM EDT
 */
#include "tools.h"
#include "classical.h"
#include "domain.h"
#include "ade.h"
#include "h5io.h"
#include "excitation.h"
#include "probes.h"
#include "step.h"
#include "mur.h"
#include "pml.h"
#include "cpml.h"
#include "dft.h"

/* mode */
unsigned int mode;
unsigned int pml_type;

/* constants */
double MU0;
double C0;
double EPSILON0;

/* discretization */
double d_x;
double d_y;
double d_z;
double d_t;

/* simulation time */
unsigned int total_timesteps;

/* grid size */
unsigned int abc_size;
unsigned int total_x;
unsigned int total_y; 
unsigned int total_z; 
unsigned int main_length_x;
unsigned int main_length_y;
unsigned int main_length_z;

/* the fields */
double ***ex;
double ***ey;
double ***ez;
double ***dx;
double ***dy;
double ***dz;
double ***bx;
double ***by;
double ***bz;
double ***hx;
double ***hy;
double ***hz;

/* the exciation sources */
double ***dipole_ex;
double ***dipole_ey;
double ***dipole_ez;
double ***dipole_hx;
double ***dipole_hy;
double ***dipole_hz;

/* domain partition */
DomainData  partition_data[7];

static void compute (int time_index)
{
    get_h();
    apply_hhards(time_index);
    get_e();
    update_mur();
    apply_ehards(time_index);
    apply_pmc();
}

int main (int argc, char** argv) 
{
    char *input_file_name;
    char *output_file_name;
    int time_index;
    int status;
    int verbosity;

    check(argc > 1, "less than two arguments; exit !");
    input_file_name = argv[1];
    output_file_name = argv[2];

    status = h5_set_file(output_file_name);
    check(status, "fail to create h5 output files");

    status = h5_get_attr(input_file_name, "settings", "number_of_timesteps", &total_timesteps);
    check(status, "fail to get h5 attributes");
    status = h5_get_attr(input_file_name, "settings", "verbosity", &verbosity);
    check(status && verbosity >= 0, "fail to get h5 attributes");
    status = h5_get_attr(input_file_name, "boundaries", "pml_type", &pml_type);
    check(status, "fail to get h5 attributes");

    status = setup_domain(input_file_name);
    check(status, "fail to setup domain");

    
    status = setup_fields(input_file_name);
    check(status, "fail to setup fields");

    status = setup_mur();
    check(status, "fail to setup abc_mur");

    status = setup_classical(input_file_name);
    check(status, "fail to setup classical mode");
    if (verbosity > 0) printf("setup classical ... DONE\n");

    if (classical != 1)
    {
        status = setup_ade(input_file_name);
        check(status, "fail to setup ade");
    }

    if (pml_type != 0) 
    {
        status = setup_cpml(input_file_name);
        check(status, "fail to setup cpml");

        if (verbosity > 0) printf("setup cpml ... DONE\n");
    }
    else  
    {
        status = setup_pml(input_file_name);
        check(status, "fail to setup pml");
        if (verbosity > 0) printf("setup pml ... DONE\n");
    }

    status = setup_input_ports(input_file_name);
    check(status, "fail to setup input ports");
    if (verbosity > 0) printf("setup excitation ... DONE\n");

    status = setup_point_sources(input_file_name);
    check(status, "fail to setup point sources");
    if (verbosity > 0) printf("setup point sources ... DONE\n");

    status = setup_hards(input_file_name);
    check(status, "fail to setup hards");
    if (verbosity > 0) printf("setup hard sources ... DONE \n");

    status = setup_output_ports(input_file_name);
    check(status, "fail to setup output ports");
    if (verbosity > 0) printf("setup probes ... DONE\n");


    status = setup_planes(input_file_name, output_file_name);
    check(status, "fail to setup planes");
    if (verbosity > 0) printf("setup animations ... DONE\n");

    status = setup_dft(input_file_name);
    check(status, "fail to setup dft");
    if (verbosity > 0) printf("setup DFT ... DONE\n");

    for (time_index = 0; time_index < total_timesteps; time_index++)
    {
        if (verbosity > 0 && time_index%verbosity == 0) printf("Time Index %06d:\n", time_index);

        excite(time_index);
        compute(time_index);
        update_ports(time_index);
        update_dft(time_index);
        update_planes(output_file_name, time_index);
    }

    write_ports(output_file_name);
    write_dft(output_file_name);
    exit(EXIT_SUCCESS);
}

