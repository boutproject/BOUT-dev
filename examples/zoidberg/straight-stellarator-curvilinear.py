#!/usr/bin/env python

# Create grids for a straight stellarator, based on 
# curvilinear grids

import numpy as np

from zoidberg import rzline
from zoidberg import field
from zoidberg import fieldtracer 
from zoidberg import poloidal_grid

#############################################################################
# Define the magnetic field
    
# Length in y after which the coils return to their starting (R,Z) locations
yperiod = 10.

magnetic_field = field.StraightStellarator(I_coil=0.3, radius = 1.0, yperiod = yperiod)

#############################################################################
# Create the inner flux surface, starting at a point at phi=0
# To do this we need to define the y locations of the poloidal points
# where we will construct grids

start_r = 0.2
start_z = 0.0

nslices = 8  # Number of poloidal slices
ycoords = np.linspace(0, yperiod, nslices)
npoints = 20  # Points per poloidal slice

# Create a field line tracer
tracer = fieldtracer.FieldTracer(magnetic_field)

# Extend the y coordinates so the tracer loops npoints times around yperiod
ycoords_all = ycoords
for i in range(1,npoints):
    ycoords_all = np.append(ycoords_all, ycoords + i*yperiod)
    
coord = tracer.follow_field_lines(start_r, start_z, ycoords_all, rtol=1e-12)
    
inner_lines = []
for i in range(nslices):
    r = coord[i::nslices,0]
    z = coord[i::nslices,1]
    line = rzline.line_from_points(r,z)
    # Re-map the points so they're approximately uniform in distance along the surface
    # Note that this results in some motion of the line
    line = line.orderByDistance()
    inner_lines.append(line)
    
# Now have a list of y coordinates (ycoords) and inner lines (inner_lines)

#############################################################################
# Generate a fixed circle for the outer boundary

outer_line = rzline.circle(R0=0.0, r=0.8)

#############################################################################
# Now have inner and outer boundaries for each poloidal slice
# Generate a grid on each poloidal slice using the elliptic grid generator

nx = 20
ny = 20

pol_slices = [ poloidal_grid.grid_elliptic(inner_line, outer_line, nx,ny, show=True) for inner_line in inner_lines ]
