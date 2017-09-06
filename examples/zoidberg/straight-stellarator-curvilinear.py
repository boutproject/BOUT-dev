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

yslice = 0
pol, ycoord = grid.getPoloidalGrid(yslice)
pol_next, ycoord_next = grid.getPoloidalGrid(yslice+1)

# Plot the points on yslice+1 as 'bx'
# Note: ravel used here only so multiple labels are not created
plt.plot(np.ravel(pol_next.R), np.ravel(pol_next.Z), 'bx', label="Grid points on slice {0}".format(yslice+1))

# Plot the forward map from yslice to yslice+1 as red 'o'
forward_R = maps['forward_R'][:,yslice,:]
forward_Z = maps['forward_Z'][:,yslice,:]
plt.plot(np.ravel(forward_R), np.ravel(forward_Z), 'ro', label="Forward map from slice {0}".format(yslice))

# Mark the points which hit the inner boundary
# These are marked with a negative x index
in_boundary = maps['forward_xt_prime'][:,yslice,:] < 0.0
plt.plot( np.ravel(forward_R[ in_boundary ]), np.ravel(forward_Z[ in_boundary ]), 'ko', label="Inner boundary points")

# Outer boundary marked with x index nx
out_boundary = maps['forward_xt_prime'][:,yslice,:] > nx-0.5
plt.plot( np.ravel(forward_R[ out_boundary ]), np.ravel(forward_Z[ out_boundary ]), 'bo', label="Outer boundary points")


plt.legend()

plt.show()
