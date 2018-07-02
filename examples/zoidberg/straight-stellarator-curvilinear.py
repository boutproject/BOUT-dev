#!/usr/bin/env python

# Create grids for a straight stellarator, based on curvilinear grids

import numpy as np

import zoidberg


#############################################################################
# Define the magnetic field
    
# Length in y after which the coils return to their starting (R,Z) locations
yperiod = 10.

magnetic_field = zoidberg.field.StraightStellarator(I_coil=0.3, radius = 1.0, yperiod = yperiod)

#############################################################################
# Create the inner flux surface, starting at a point at phi=0
# To do this we need to define the y locations of the poloidal points
# where we will construct grids

start_r = 0.2
start_z = 0.0

nslices = 8  # Number of poloidal slices
ycoords = np.linspace(0, yperiod, nslices)
npoints = 20  # Points per poloidal slice

rzcoord, ycoords = zoidberg.fieldtracer.trace_poincare(magnetic_field, start_r, start_z, yperiod,
                                                     y_slices=ycoords, revs=npoints)
    
inner_lines = []
for i in range(nslices):
    r = rzcoord[:,i,0]
    z = rzcoord[:,i,1]
    line = zoidberg.rzline.line_from_points(r,z)
    # Re-map the points so they're approximately uniform in distance along the surface
    # Note that this results in some motion of the line
    line = line.equallySpaced()
    inner_lines.append(line)
    
# Now have a list of y coordinates (ycoords) and inner lines (inner_lines)

#############################################################################
# Generate a fixed circle for the outer boundary

outer_line = zoidberg.rzline.circle(R0=0.0, r=0.8)

#############################################################################
# Now have inner and outer boundaries for each poloidal grid
# Generate a grid on each poloidal slice using the elliptic grid generator

nx = 20
ny = 20

pol_grids = [ zoidberg.poloidal_grid.grid_elliptic(inner_line, outer_line, nx,ny) for inner_line in inner_lines ]

#############################################################################
# Create a grid, then calculate forward and backward maps

grid = zoidberg.grid.Grid( pol_grids, ycoords, yperiod, yperiodic=True)

maps = zoidberg.make_maps(grid, magnetic_field)

#############################################################################
# Write grid file

filename = "stellarator.fci.nc"

print("Writing to grid file '{0}'".format(filename))
zoidberg.write_maps(grid, magnetic_field, maps, gridfile=filename, new_names=False, metric2d=True)

#############################################################################
# Plot maps

zoidberg.plot.plot_forward_map(grid, maps)
