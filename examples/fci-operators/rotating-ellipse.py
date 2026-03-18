#!/usr/bin/env python

# Create grids for a rotating ellipse stellarator, based on curvilinear grids

import numpy as np

import zoidberg

#############################################################################
# Define the magnetic field

# Length in y after which the coils return to their starting (R,Z) locations
yperiod = 2.0 * np.pi / 5

xcentre = 2.0

magnetic_field = zoidberg.field.RotatingEllipse(
    xcentre=xcentre,  # Major radius of axis [m]
    radius=0.8,  # Minor radius of the coils [m]
    I_coil=0.4,  # Coil current
    yperiod=yperiod,
    Btor=1.0,  # Toroidal magnetic field
)

#############################################################################
# Create the inner flux surface, starting at a point at phi=0
# To do this we need to define the y locations of the poloidal points
# where we will construct grids

start_r = xcentre + 0.4
start_z = 0.0

nslices = 8  # Number of poloidal slices
ycoords = np.linspace(0, yperiod, nslices)
npoints = 20  # Points per poloidal slice

rzcoord, ycoords = zoidberg.fieldtracer.trace_poincare(
    magnetic_field, start_r, start_z, yperiod, y_slices=ycoords, revs=npoints, nover=1
)

inner_lines = []
for i in range(nslices):
    r = rzcoord[:, i, 0, 0]
    z = rzcoord[:, i, 0, 1]
    line = zoidberg.rzline.line_from_points(r, z)
    # Re-map the points so they're approximately uniform in distance along the surface
    # Note that this results in some motion of the line
    line = line.equallySpaced()
    inner_lines.append(line)

# Now have a list of y coordinates (ycoords) and inner lines (inner_lines)

#############################################################################
# Generate a fixed circle for the outer boundary

outer_line = zoidberg.rzline.circle(R0=xcentre, r=0.6)

#############################################################################
# Now have inner and outer boundaries for each poloidal grid
# Generate a grid on each poloidal slice using the elliptic grid generator

nx = 20
nz = 20

pol_grids = [
    zoidberg.poloidal_grid.grid_elliptic(inner_line, outer_line, nx, nz)
    for inner_line in inner_lines
]

#############################################################################
# Create a grid, then calculate forward and backward maps

grid = zoidberg.grid.Grid(pol_grids, ycoords, yperiod, yperiodic=True)

maps = zoidberg.make_maps(grid, magnetic_field)

#############################################################################
# Interpolation weights

weight_vars = zoidberg.weights.calculate_parallel_map_operators(maps)
maps.update(weight_vars)

#############################################################################
# Write grid file

filename = f"rotating-ellipse-{nx}x{nslices}x{nz}.fci.nc"

print(f"Writing to grid file '{filename}'")
zoidberg.write_maps(
    grid, magnetic_field, maps, gridfile=filename, new_names=False, metric2d=False
)
