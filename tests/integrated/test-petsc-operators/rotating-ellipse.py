#!/usr/bin/env python

import argparse
import numpy as np
import zoidberg

parser = argparse.ArgumentParser(
    prog="RotatingEllipse",
    description="Create grids for a rotating ellipse stellarator, based on curvilinear grids.",
)
parser.add_argument(
    "-nx", type=int, default=20, help="Number of radial (x) cells including boundaries"
)
parser.add_argument(
    "-ny", "--nslices", type=int, default=8, help="Number of parallel (y) slices"
)
parser.add_argument("-nz", type=int, default=20, help="Number of poloidal (z) cells")
parser.add_argument(
    "-ns",
    "--nsamples-per-dim",
    type=int,
    default=1,
    help="Number of samples per dimension",
)
parser.add_argument("-pb", "--partial-boundaries", action="store_true", default=True)
parser.add_argument("-o", "--output", type=str, default="")
args = parser.parse_args()

if args.output == "":
    args.output = f"rotating-ellipse-{args.nx}x{args.nslices}x{args.nz}-s{args.nsamples_per_dim}{'-npb' if not args.partial_boundaries else ''}.fci.nc"

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

ycoords = np.linspace(0, yperiod, args.nslices, endpoint=False)
npoints = 20  # Points per poloidal slice

rzcoord, ycoords = zoidberg.fieldtracer.trace_poincare(
    magnetic_field, start_r, start_z, yperiod, y_slices=ycoords, revs=npoints, nover=1
)

inner_lines = []
for i in range(args.nslices):
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

pol_grids = [
    zoidberg.poloidal_grid.grid_elliptic(inner_line, outer_line, args.nx, args.nz)
    for inner_line in inner_lines
]

#############################################################################
# Create a grid, then calculate forward and backward maps

grid = zoidberg.grid.Grid(pol_grids, ycoords, yperiod, yperiodic=True)

maps = zoidberg.make_maps(grid, magnetic_field, samples_per_dim=args.nsamples_per_dim)

maps = zoidberg.weights.modify_maps(maps, partial_boundaries=args.partial_boundaries)

#############################################################################
# Write grid file

print(f"Writing to grid file '{args.output}'")
zoidberg.write_maps(
    grid, magnetic_field, maps, gridfile=args.output, new_names=False, metric2d=False
)
