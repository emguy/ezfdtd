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
 *     File Name : matrix.c
 * Last Modified : Tue 24 Sep 2013 05:58:55 PM EDT
 */
#include "tools.h"
#include "mem.h"

void ***mem3(int type, size_t xmax, size_t ymax, size_t zmax) 
{
    unsigned int  x, y;
    void ***ppp;

    ppp = NULL;

    switch (type)
    {
        case type_int:
            ppp = malloc(xmax * sizeof(int **));
            inspect(ppp, KCYN "malloc failed for %d bytes", xmax * sizeof(int **));

            for (x = 0; x < xmax; x++) 
            {
                ppp[x] = malloc(ymax * sizeof(int *));
                inspect(ppp[x], KCYN "malloc failed for %d bytes", ymax * sizeof(int *));

                for (y = 0; y < ymax; y++) 
                {
                    ppp[x][y] = (int *)calloc(zmax, sizeof(int));
                    inspect(ppp[x][y], KCYN "malloc failed for %d bytes", zmax * sizeof(int));
                }
            }
            break;

        case type_double:
            ppp = malloc(xmax * sizeof(double **));
            inspect(ppp, KCYN "malloc failed for %d bytes", xmax * sizeof(double **));
            for (x = 0; x < xmax; x++) 
            {
                ppp[x] = malloc(ymax * sizeof(double *));
                inspect(ppp[x], KCYN "malloc failed for %d bytes", ymax * sizeof(double *));

                for (y = 0; y < ymax; y++) 
                {
                    ppp[x][y] = (double *)calloc(zmax, sizeof(double));
                    inspect(ppp[x][y], KCYN "malloc failed for %d bytes", zmax * sizeof(double));
                }
            }
            break;
    }

    return(ppp);
}

void** mem2(int type, size_t xmax, size_t ymax) 
{
    unsigned int  x;
    void  **pp;
    pp = NULL;

    switch (type)
    {
        case type_int:
            pp = malloc(xmax * sizeof(int *));
            inspect(pp, KCYN "malloc failed for %d bytes", xmax * sizeof(int *));

            for (x = 0; x < xmax; x++) 
            {
                pp[x] = (int *)calloc(ymax, sizeof(int));
                inspect(pp[x], KCYN "malloc failed for %d bytes", ymax * sizeof(int));
            }
            break;

        case type_double:
            pp = malloc(xmax * sizeof(double **));
            inspect(pp, KCYN "malloc failed for %d bytes", xmax * sizeof(double **));
            for (x = 0; x < xmax; x++) 
            {
                pp[x] = (double *)calloc(ymax, sizeof(double));
                inspect(pp[x], KCYN "malloc failed for %d bytes", ymax * sizeof(double));
            }
            break;
    }

    return(pp);
}

void *mem1(int type, size_t xmax)
{
    void  *p;
    p = NULL;

    switch (type)
    {
        case type_int:
            p = (int *)calloc(xmax, sizeof(int));
            inspect(p, KCYN "malloc failed for %d bytes", xmax * sizeof(int));
            break;

        case type_double:
            p = (double *)calloc(xmax, sizeof(double));
            inspect(p, KCYN "malloc failed for %d bytes", xmax * sizeof(double));
            break;
    }

    return(p);
}

void **contiguous_mem2 (int type, size_t xmax, size_t ymax) 
{
    unsigned int x, xy;
    void **pp, *p;
    pp = NULL;

    xy = xmax * ymax;

    switch (type)
    {
        case type_int:
            p = (int *)calloc(xy, sizeof(int));
            inspect(p, KCYN "malloc failed for %d bytes", xy * sizeof(int));
            pp = malloc(xmax * sizeof(int *));
            inspect(pp, KCYN "malloc failed for %d bytes", xmax * sizeof(int *));
            for (x = 0; x < xmax; x++) pp[x] = &((int *)p)[ymax * x];
            break;

        case type_double:
            p = (double *)calloc(xy, sizeof(double));
            inspect(p, KCYN "malloc failed for %d bytes", xy * sizeof(double));
            pp = malloc(xmax * sizeof(double *));
            inspect(pp, KCYN "malloc failed for %d bytes", xmax * sizeof(double *));
            for (x = 0; x < xmax; x++) pp[x] = &((double *)p)[ymax * x];
            break;
    }
    return(pp);
}

void ***contiguous_mem3 (int type, size_t xmax, size_t ymax, size_t zmax) 
{
    unsigned int x, y, xy, xyz;
    void ***ppp, **pp, *p;

    ppp = NULL;

    xyz = xmax * ymax * zmax;
    xy = xmax * ymax;

    switch (type)
    {
        case type_int:
            p = (int *)calloc(xyz, sizeof(int));
            inspect(p, KCYN "malloc failed for %d bytes", xyz * sizeof(int));
            pp = malloc(xy * sizeof(int *));
            inspect(pp, KCYN "malloc failed for %d bytes", xy * sizeof(int *));
            ppp = malloc(xmax * sizeof(int **));
            inspect(ppp, KCYN "malloc failed for %d bytes", xmax * sizeof(int **));
            break;
            for (x = 0; x < xmax; x++)
            {
                ppp[x] = pp + ymax * x;
                for (y = 0; y < ymax; y++)
                    ppp[x][y] = (int *)p + ymax * zmax * x + zmax * y;
            }
            break;

        case type_double:
            p = (double *)calloc(xyz, sizeof(double));
            inspect(p, KCYN "malloc failed for %d bytes", xyz * sizeof(double));
            pp = malloc(xy * sizeof(double *));
            inspect(pp, KCYN "malloc failed for %d bytes", xy * sizeof(double *));
            ppp = malloc(xmax * sizeof(double **));
            inspect(ppp, KCYN "malloc failed for %d bytes", xmax * sizeof(double **));
            break;
            for (x = 0; x < xmax; x++)
            {
                ppp[x] = pp + ymax * x;
                for (y = 0; y < ymax; y++)
                    ppp[x][y] = (int *)p + ymax * zmax * x + zmax * y;
            }
    }
    return(ppp);
}

void free_mem3 (void ***mem, int xmax, int ymax)
{
    unsigned int x, y;
    for (x = 0; x < xmax; x++)
    {
        for (y = 0; y < ymax; y++)
            free(mem[x][y]);
        free(mem[x]);
    }
    free(mem);
}

void free_mem2 (void **mem, int xmax)
{
    unsigned int x;
    for (x = 0; x < xmax; x++)
        free(mem[x]);
    free(mem);
}

void free_mem1 (void *mem)
{
    free(mem);
}

void free_contiguous_mem2 (void **mem)
{
    free(mem[0]);
    free(mem);
}

void free_contiguous_mem3 (void ***mem)
{
    free(mem[0][0]);
    free(mem[0]);
    free(mem);
}
