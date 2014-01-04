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
 *     File Name : h5io.c
 * Last Modified : Wed 25 Sep 2013 03:55:29 PM EDT
 */

#include "tools.h"
#include "mem.h"
#include <hdf5.h>
#include "h5io.h"

#define H5_DEBUG 1

void ***h5_load3(const char *file_name, const char *dset_name, size_t xmax, size_t ymax, size_t zmax)
{
    hid_t h5_file, h5_space, h5_dset, h5_type;
    hsize_t dims[3];
    herr_t status;
    void *p;
    void ***ppp;
    int type;
    int x, y, z;

    type = -1;
    ppp = NULL;

    h5_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);
    h5_dset = H5Dopen(h5_file, dset_name, H5P_DEFAULT);
    inspect(h5_dset > 0, "fail to access dset \"%s\" (%s)", dset_name, file_name);
    h5_space = H5Dget_space(h5_dset);
    inspect(H5Sget_simple_extent_ndims(h5_space) == 3, "dimensionality mismatch (%d, %d) in \"%s\" (%s)", H5Sget_simple_extent_ndims(h5_space), 3, dset_name, file_name);
    H5Sget_simple_extent_dims(h5_space, dims, NULL);
    inspect(dims[0]==xmax || dims[1]==ymax || dims[2]==zmax, "dimension mismatch (%dx%dx%d, %dx%dx%d) in \"%s\" (%s)", dims[0], dims[1], dims[2], xmax, ymax, zmax, dset_name, file_name);
    h5_type = H5Tget_native_type(H5Dget_type(h5_dset), H5T_DIR_ASCEND);
    inspect(h5_type > 0, "unkown datatype in \"%s\" from %s", dset_name, file_name);

    if (H5Tequal(h5_type, H5T_NATIVE_INT) || H5Tequal(h5_type, H5T_NATIVE_LONG) || H5Tequal(h5_type, H5T_NATIVE_LLONG) || H5Tequal(h5_type, H5T_NATIVE_SHORT))
        type = type_int;
    else if (H5Tequal(h5_type, H5T_NATIVE_DOUBLE) || H5Tequal(h5_type, H5T_NATIVE_FLOAT) || H5Tequal(h5_type, H5T_NATIVE_LDOUBLE))
        type = type_double;
    inspect(type > 0, "invalid datatype in \"%s\" (%s)", dset_name, file_name);

    p = mem1(type, xmax * ymax * zmax);
    inspect(p, "fail to allocate memory");
    ppp = mem3(type, xmax, ymax, zmax);
    inspect(ppp, "fail to allocate memory");

    switch (type)
    {
        case type_int:
            status = H5Dread(h5_dset, H5T_NATIVE_INT, h5_space, H5S_ALL, H5P_DEFAULT, p);
            break;
        case type_double:
            status = H5Dread(h5_dset, H5T_NATIVE_DOUBLE, h5_space, H5S_ALL, H5P_DEFAULT, p);
            break;
        default:
            status = -1;
            break;
    }
    inspect(status >= 0," fail to read data from \"%s\" (%s)", dset_name, file_name);

    H5Tclose(h5_type);
    H5Sclose(h5_space);
    H5Dclose(h5_dset);
    H5Fclose(h5_file);

    for (x = 0; x < xmax; x++)
        for (y = 0; y < ymax; y++)
            for (z = 0; z < zmax; z++)
                switch (type)
                {
                    case type_int:
                        ((int ***)ppp)[x][y][z] = ((int *)p)[x*ymax*zmax + y*zmax + z];
                        break;
                    case type_double:
                        ((double ***)ppp)[x][y][z] = ((double *)p)[x*ymax*zmax + y*zmax + z];
                }
    free(p);
    return ppp;
}

void **h5_load2(const char *file_name, const char *dset_name, size_t xmax, size_t ymax)
{
    hid_t h5_file, h5_space, h5_dset, h5_type;
    hsize_t dims[2];
    herr_t status;
    void *p;
    void **pp;
    int type;
    int x, y;

    type = -1;
    pp = NULL;

    h5_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);
    h5_dset = H5Dopen(h5_file, dset_name, H5P_DEFAULT);
    inspect(h5_dset > 0, "fail to access dset \"%s\" (%s)", dset_name, file_name);
    h5_space = H5Dget_space(h5_dset);
    inspect(H5Sget_simple_extent_ndims(h5_space) == 2, "dimensionality mismatch (%d, %d) in \"%s\" (%s)", H5Sget_simple_extent_ndims(h5_space), 2, dset_name, file_name);
    H5Sget_simple_extent_dims(h5_space, dims, NULL);
    inspect(dims[0]==xmax || dims[1]==ymax, "dimension mismatch (%dx%d, %dx%d) in \"%s\" (%s)", dims[0], dims[1], xmax, ymax, dset_name, file_name);
    h5_type = H5Tget_native_type(H5Dget_type(h5_dset), H5T_DIR_ASCEND);
    inspect(h5_type > 0, "unkown datatype in \"%s\" from %s", dset_name, file_name);

    if (H5Tequal(h5_type, H5T_NATIVE_INT) || H5Tequal(h5_type, H5T_NATIVE_LONG) || H5Tequal(h5_type, H5T_NATIVE_LLONG) || H5Tequal(h5_type, H5T_NATIVE_SHORT))
        type = type_int;
    else if (H5Tequal(h5_type, H5T_NATIVE_DOUBLE) || H5Tequal(h5_type, H5T_NATIVE_FLOAT) || H5Tequal(h5_type, H5T_NATIVE_LDOUBLE))
        type = type_double;
    inspect(type > 0, "invalid datatype in \"%s\" (%s)", dset_name, file_name);

    p = mem1(type, xmax * ymax);
    inspect(p, "fail to allocate memory");
    pp = mem2(type, xmax, ymax);
    inspect(pp, "fail to allocate memory");

    switch (type)
    {
        case type_int:
            status = H5Dread(h5_dset, H5T_NATIVE_INT, h5_space, H5S_ALL, H5P_DEFAULT, p);
            break;
        case type_double:
            status = H5Dread(h5_dset, H5T_NATIVE_DOUBLE, h5_space, H5S_ALL, H5P_DEFAULT, p);
            break;
        default:
            status = -1;
            break;
    }
    inspect(status >= 0," fail to read data from \"%s\" (%s)", dset_name, file_name);

    H5Tclose(h5_type);
    H5Sclose(h5_space);
    H5Dclose(h5_dset);
    H5Fclose(h5_file);

    for (x = 0; x < xmax; x++)
        for (y = 0; y < ymax; y++)
            switch (type)
            {
                case type_int:
                    ((int **)pp)[x][y] = ((int *)p)[x*ymax + y];
                    break;
                case type_double:
                    ((double **)pp)[x][y] = ((double *)p)[x*ymax + y];
            }

    free(p);
    return pp;
}

void **h5_load1(const char *file_name, const char *dset_name, size_t xmax)
{
    hid_t h5_file, h5_space, h5_dset, h5_type;
    hsize_t dims[1];
    herr_t status;
    void *p;
    int type;

    status = -1;
    type = -1;
    p = NULL;

    h5_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);
    h5_dset = H5Dopen(h5_file, dset_name, H5P_DEFAULT);
    inspect(h5_dset > 0, "fail to access dset \"%s\" (%s)", dset_name, file_name);
    h5_space = H5Dget_space(h5_dset);
    inspect(H5Sget_simple_extent_ndims(h5_space) == 1, "dimensionality mismatch (%d, %d) in \"%s\" (%s)", H5Sget_simple_extent_ndims(h5_space), 1, dset_name, file_name);
    H5Sget_simple_extent_dims(h5_space, dims, NULL);
    inspect(dims[0]==xmax, "dimension mismatch in \"%s\" (%s)", dset_name, file_name);
    inspect(dims[0]==xmax, "dimension mismatch (%dx%d, %dx%d) in \"%s\" (%s)", dims[0], xmax, dset_name, file_name);
    h5_type = H5Tget_native_type(H5Dget_type(h5_dset), H5T_DIR_ASCEND);
    inspect(h5_type > 0, "unkown datatype in \"%s\" from %s", dset_name, file_name);

    if (H5Tequal(h5_type, H5T_NATIVE_INT) || H5Tequal(h5_type, H5T_NATIVE_LONG) || H5Tequal(h5_type, H5T_NATIVE_LLONG) || H5Tequal(h5_type, H5T_NATIVE_SHORT))
        type = type_int;
    else if (H5Tequal(h5_type, H5T_NATIVE_DOUBLE) || H5Tequal(h5_type, H5T_NATIVE_FLOAT) || H5Tequal(h5_type, H5T_NATIVE_LDOUBLE))
        type = type_double;
    inspect(type > 0, "invalid datatype in \"%s\" (%s)", dset_name, file_name);

    p = mem1(type, xmax);
    inspect(p, "fail to allocate memory");

    switch (type)
    {
        case type_int:
            status = H5Dread(h5_dset, H5T_NATIVE_INT, h5_space, H5S_ALL, H5P_DEFAULT, p);
            break;
        case type_double:
            status = H5Dread(h5_dset, H5T_NATIVE_DOUBLE, h5_space, H5S_ALL, H5P_DEFAULT, p);
            break;
        default:
            status = -1;
            break;
    }
    inspect(status >= 0,"fail to read data from \"%s\" (%s)", dset_name, file_name);

    H5Tclose(h5_type);
    H5Sclose(h5_space);
    H5Dclose(h5_dset);
    H5Fclose(h5_file);

    return p;
}

int h5_get_attr(const char *file_name, const char *grp_name, const char *attr_name, void *value)
{
    hid_t h5_file, h5_grp, h5_attr, h5_type;
    herr_t status;
    int type;

    type = -1;
    h5_type = -1;

    h5_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);
    h5_grp = H5Gopen2(h5_file, grp_name, H5P_DEFAULT);
    inspect(h5_grp >  0, "fail to access node \"%s\" (%s)", grp_name, file_name);
    h5_attr = H5Aopen_by_name(h5_file, grp_name, attr_name, H5P_DEFAULT, H5P_DEFAULT);
    inspect(h5_attr > 0, "fail to access attribute \"%s.%s\" (%s)",grp_name, attr_name, file_name);
    h5_type = H5Tget_native_type(H5Aget_type(h5_attr), H5T_DIR_ASCEND);
    inspect(h5_type > 0, "unknown datatype in attr \"%s.%s\" (%s)", grp_name, attr_name, file_name);

    /*

    if(H5_DEBUG == 1)
    {
        printf("H5T_NATIVE_CHAR    :%d\n", H5Tequal(h5_type, H5T_NATIVE_CHAR   )   );
        printf("H5T_NATIVE_SHORT   :%d\n", H5Tequal(h5_type, H5T_NATIVE_SHORT  )   );
        printf("H5T_NATIVE_INT     :%d\n", H5Tequal(h5_type, H5T_NATIVE_INT    )   ); 
        printf("H5T_NATIVE_LONG    :%d\n", H5Tequal(h5_type, H5T_NATIVE_LONG   )   );
        printf("H5T_NATIVE_LLONG   :%d\n", H5Tequal(h5_type, H5T_NATIVE_LLONG  )   );
        printf("H5T_NATIVE_UCHAR   :%d\n", H5Tequal(h5_type, H5T_NATIVE_UCHAR  )   );
        printf("H5T_NATIVE_USHORT  :%d\n", H5Tequal(h5_type, H5T_NATIVE_USHORT )   );
        printf("H5T_NATIVE_UINT    :%d\n", H5Tequal(h5_type, H5T_NATIVE_UINT   )   );
        printf("H5T_NATIVE_ULONG   :%d\n", H5Tequal(h5_type, H5T_NATIVE_ULONG  )   );
        printf("H5T_NATIVE_ULLONG  :%d\n", H5Tequal(h5_type, H5T_NATIVE_ULLONG )   );
        printf("H5T_NATIVE_FLOAT   :%d\n", H5Tequal(h5_type, H5T_NATIVE_FLOAT  )   );
        printf("H5T_NATIVE_DOUBLE  :%d\n", H5Tequal(h5_type, H5T_NATIVE_DOUBLE )   );
        printf("H5T_NATIVE_LDOUBLE :%d\n", H5Tequal(h5_type, H5T_NATIVE_LDOUBLE)   );
        printf("H5T_NATIVE_B8      :%d\n", H5Tequal(h5_type, H5T_NATIVE_B8     )   );
        printf("H5T_NATIVE_B16     :%d\n", H5Tequal(h5_type, H5T_NATIVE_B16    )   );
        printf("H5T_NATIVE_B32     :%d\n", H5Tequal(h5_type, H5T_NATIVE_B32    )   );
        printf("H5T_NATIVE_B64     :%d\n", H5Tequal(h5_type, H5T_NATIVE_B64    )   );
    }
    */


    if (H5Tequal(h5_type, H5T_NATIVE_INT) || H5Tequal(h5_type, H5T_NATIVE_LONG) || H5Tequal(h5_type, H5T_NATIVE_LLONG) || H5Tequal(h5_type, H5T_NATIVE_SHORT))
        type = type_int;
    else if (H5Tequal(h5_type, H5T_NATIVE_DOUBLE) || H5Tequal(h5_type, H5T_NATIVE_FLOAT) || H5Tequal(h5_type, H5T_NATIVE_LDOUBLE))
        type = type_double;
    inspect(type > 0, "invalid datatype in attr \"%s.%s\" (%s)", grp_name, attr_name, file_name);

    switch (type)
    {
        case type_int:
            status = H5Aread(h5_attr, H5T_NATIVE_INT, value);
            break;
        case type_double:
            status = H5Aread(h5_attr, H5T_NATIVE_DOUBLE, value);
            break;
        default:
            status = -1;
            break;
    }
    inspect(status >= 0,"fail to read attribute \"%s.%s\" (%s)", grp_name, attr_name, file_name);

    H5Tclose(h5_type);
    H5Aclose(h5_attr);
    H5Gclose(h5_grp);
    H5Fclose(h5_file);

    return 1;
}

int h5_set_file(const char *file_name)
{
    hid_t           h5_file;
    herr_t          status;
    htri_t          avail;
    unsigned int    filter_info;

    avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);
    inspect(avail, "gzip filter not available");

    status = H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
    inspect(status >= 0, "gzip info not available");
    inspect( (filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) &&
            (filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED), 
            "gzip filter not available for encoding and decoding");

    h5_file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to create file %s", file_name);
    H5Fclose(h5_file);

    return 1;
}

int h5_write2(const char *file_name, const char *dset_name, int type, void **data, size_t xmax, size_t ymax)
{
    hid_t           h5_file, h5_space, dcpl;
    hid_t           h5_dset;
    herr_t          status;
    hsize_t         dims[2] = {xmax, ymax};

    void *p;
    int x, y;

    h5_file = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    h5_space = H5Screate_simple(2, dims, NULL);

    p = mem1(type, xmax * ymax);
    inspect(p, "fail to allocate memory");

    for (x = 0; x < xmax; x++)
        for (y = 0; y < ymax; y++)
            switch (type)
            {
                case type_int:
                    ((int *)p)[x*ymax + y] = ((int **)data)[x][y];
                    break;
                case type_double:
                    ((double *)p)[x*ymax + y] = ((double **)data)[x][y];
                    break;
            }

    switch (type)
    {
        case type_int:
            h5_dset = H5Dcreate(h5_file, dset_name, H5T_NATIVE_INT, h5_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            status = H5Dwrite(h5_dset, H5T_NATIVE_INT, h5_space, h5_space, H5P_DEFAULT, p);
            break;
        case type_double:
            h5_dset = H5Dcreate(h5_file, dset_name, H5T_NATIVE_DOUBLE, h5_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            status = H5Dwrite(h5_dset, H5T_NATIVE_DOUBLE, h5_space, h5_space, H5P_DEFAULT, p);
            break;
        default:
            status = -1;
            break;
    }
    inspect(status >= 0, "fail to write dataset \"%s\" to file (%s)", dset_name, file_name);

    H5Pclose(dcpl);
    H5Dclose(h5_dset);
    H5Sclose(h5_space);
    H5Fclose(h5_file);

    return 1;
}

int h5_empty3(const char *file_name, const char *dset_name, int type, size_t xmax, size_t ymax, size_t zmax)
{
    hsize_t         dims[3] = {xmax, ymax, zmax};
    hsize_t         chunk[3] = {xmax, ymax, 1};
    hid_t           h5_file, h5_space, h5_dset, dcpl;

    h5_file = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);

    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    h5_space = H5Screate_simple(3, dims, NULL);
    H5Pset_deflate(dcpl,9);
    H5Pset_chunk(dcpl, 3, chunk);

    switch (type)
    {
        case type_int:
            h5_dset = H5Dcreate(h5_file, dset_name, H5T_NATIVE_INT, h5_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            break;
        case type_double:
            h5_dset = H5Dcreate(h5_file, dset_name, H5T_NATIVE_DOUBLE, h5_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            break;
        default:
            h5_dset = -1;
            break;
    }
    inspect(h5_dset >= 0, "fail to create dataset \"%s\" to file (%s)", dset_name, file_name);

    H5Pclose(dcpl);
    H5Dclose(h5_dset);
    H5Sclose(h5_space);
    H5Fclose(h5_file);

    return 1;
}

int h5_slab(const char *file_name, const char *dset_name, int type, void **buffer, size_t xmax, size_t ymax, int index)
{
    hid_t           h5_file, h5_space, slab_space, h5_dset;
    herr_t          status;
    hsize_t         slab_dims[2] = {xmax, ymax};
    hsize_t         count[3] = {xmax, ymax, 1};
    hsize_t         offset[3] = {0, 0, index};

    h5_file = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    inspect(h5_file > 0, "fail to access file %s", file_name);
    h5_dset = H5Dopen(h5_file, dset_name, H5P_DEFAULT);
    inspect(h5_dset > 0, "fail to access dset \"%s\" (%s)", dset_name, file_name);
    h5_space = H5Dget_space(h5_dset);
    slab_space = H5Screate_simple(2, slab_dims, NULL);

    status = H5Sselect_hyperslab(h5_space, H5S_SELECT_SET, offset, NULL, count, NULL);
    inspect(status >= 0, "fail to access the hyperslab on dset \"%s\" (%s)", dset_name, file_name);

    switch (type)
    {
        case type_int:
            status = H5Dwrite(h5_dset, H5T_NATIVE_INT, slab_space, h5_space, H5P_DEFAULT, buffer[0]);
            break;
        case type_double:
            status = H5Dwrite(h5_dset, H5T_NATIVE_DOUBLE, slab_space, h5_space, H5P_DEFAULT, buffer[0]);
            break;
        default:
            status = -1;
    }
    inspect(h5_dset >= 0, "fail to write to the dataset \"%s\" (%s)", dset_name, file_name);

    status = H5Dclose(h5_dset);
    status = H5Sclose(slab_space);
    status = H5Sclose(h5_space);
    status = H5Fclose(h5_file);

    return 1;
}

