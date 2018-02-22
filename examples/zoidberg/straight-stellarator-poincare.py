#!/usr/bin/env python

# Create grids for a straight stellarator, based on 
# curvilinear grids

import numpy as np
import matplotlib.pyplot as plt

import zoidberg

#############################################################################
# Define the magnetic field
    
# Length in y after which the coils return to their starting (R,Z) locations
yperiod = 10.

magnetic_field = zoidberg.field.StraightStellarator(I_coil=0.4, radius = 1.0, yperiod = yperiod)

start_r = 0.3
start_z = 0.0

nslices = 4  # Number of poloidal slices
ycoords = np.linspace(0, yperiod, nslices)
npoints = 60  # Points per poloidal slice

# Create a field line tracer
tracer = zoidberg.fieldtracer.FieldTracer(magnetic_field)
#tracer = zoidberg.fieldtracer.FieldTracerReversible(magnetic_field)

# Extend the y coordinates so the tracer loops npoints times around yperiod
ycoords_all = ycoords
for i in range(1,npoints):
    ycoords_all = np.append(ycoords_all, ycoords + i*yperiod)
    
coord = tracer.follow_field_lines(start_r, start_z, ycoords_all, rtol=1e-12)

for i in range(nslices):
    r = coord[i::nslices,0]
    z = coord[i::nslices,1]
    plt.plot(r,z,'o')

#i = 1
#plt.plot(coord[i::nslices,0], coord[i::nslices,1],'o')

plt.show()

    
