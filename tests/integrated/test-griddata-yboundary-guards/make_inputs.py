#! /usr/bin/env python3

import numpy
from netCDF4 import Dataset
import os.path

# make a grid with topology of a 'double null' case

nx = 4
ny = 24
blocksize = ny/6

for n_yguards in [0, 1, 2]:
    datadir = "data-doublenull-" + str(n_yguards)
    gridname = "grid-doublenull-" + str(n_yguards) + ".nc"

    with Dataset(os.path.join(datadir,gridname), 'w') as gridfile:
        gridfile.createDimension('x', nx)
        gridfile.createDimension('y', ny + 4*n_yguards)

        gridfile.createVariable('nx', numpy.int32)
        gridfile['nx'][...] = nx

        gridfile.createVariable('ny', numpy.int32)
        gridfile['ny'][...] = ny

        gridfile.createVariable('y_boundary_guards', numpy.int32)
        gridfile['y_boundary_guards'][...] = n_yguards

        gridfile.createVariable('MXG', numpy.int32)
        gridfile['MXG'][...] = 1

        gridfile.createVariable('MYG', numpy.int32)
        gridfile['MYG'][...] = 2 if n_yguards==0 else n_yguards

        gridfile.createVariable('ixseps1', numpy.int32)
        gridfile['ixseps1'][...] = nx//2 - 1

        gridfile.createVariable('ixseps2', numpy.int32)
        gridfile['ixseps2'][...] = nx//2  - 1

        gridfile.createVariable('jyseps1_1', numpy.int32)
        gridfile['jyseps1_1'][...] = blocksize - 1

        gridfile.createVariable('jyseps2_1', numpy.int32)
        gridfile['jyseps2_1'][...] = 2*blocksize - 1

        gridfile.createVariable('ny_inner', numpy.int32)
        gridfile['ny_inner'][...] = 3*blocksize

        gridfile.createVariable('jyseps1_2', numpy.int32)
        gridfile['jyseps1_2'][...] = 4*blocksize - 1

        gridfile.createVariable('jyseps2_2', numpy.int32)
        gridfile['jyseps2_2'][...] = 5*blocksize - 1

        testdata = numpy.zeros([nx, ny + 4*n_yguards])
        testdata[:,:] = numpy.arange(ny + 4*n_yguards)[numpy.newaxis,:]
        gridfile.createVariable('test', float, ('x', 'y'))
        gridfile['test'][...] = testdata

for n_yguards in [0, 1, 2]:
    datadir = "data-singlenull-" + str(n_yguards)
    gridname = "grid-singlenull-" + str(n_yguards) + ".nc"

    with Dataset(os.path.join(datadir,gridname), 'w') as gridfile:
        gridfile.createDimension('x', nx)
        gridfile.createDimension('y', ny + 2*n_yguards)

        gridfile.createVariable('nx', numpy.int32)
        gridfile['nx'][...] = nx

        gridfile.createVariable('ny', numpy.int32)
        gridfile['ny'][...] = ny

        gridfile.createVariable('y_boundary_guards', numpy.int32)
        gridfile['y_boundary_guards'][...] = n_yguards

        gridfile.createVariable('MXG', numpy.int32)
        gridfile['MXG'][...] = 1

        gridfile.createVariable('MYG', numpy.int32)
        gridfile['MYG'][...] = 2 if n_yguards==0 else n_yguards

        gridfile.createVariable('ixseps1', numpy.int32)
        gridfile['ixseps1'][...] = nx//2 - 1

        gridfile.createVariable('ixseps2', numpy.int32)
        gridfile['ixseps2'][...] = nx//2  - 1

        gridfile.createVariable('jyseps1_1', numpy.int32)
        gridfile['jyseps1_1'][...] = blocksize - 1

        gridfile.createVariable('jyseps2_1', numpy.int32)
        gridfile['jyseps2_1'][...] = ny//2

        gridfile.createVariable('ny_inner', numpy.int32)
        gridfile['ny_inner'][...] = ny//2

        gridfile.createVariable('jyseps1_2', numpy.int32)
        gridfile['jyseps1_2'][...] = ny//2

        gridfile.createVariable('jyseps2_2', numpy.int32)
        gridfile['jyseps2_2'][...] = 5*blocksize - 1

        testdata = numpy.zeros([nx, ny + 2*n_yguards])
        testdata[:,:] = numpy.arange(ny + 2*n_yguards)[numpy.newaxis,:]
        gridfile.createVariable('test', float, ('x', 'y'))
        gridfile['test'][...] = testdata

exit(0)
