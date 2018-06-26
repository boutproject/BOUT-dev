#!/usr/bin/env python

import zoidberg
import matplotlib.pyplot as plt

nx = 5
ny = 6
nz = 7

# Create magnetic field
magnetic_field = zoidberg.field.StraightStellarator(radius = 1.0)

# Create a rectangular grid in (x,y,z)
rectangle = zoidberg.grid.rectangular_grid(nx,ny,nz,
                                           Lx = 1.0, Lz = 1.0, Ly = 10.0,
                                           yperiodic = True)

# Here both the field and and grid are centred at (x,z) = (0,0)
# and the rectangular grid here fits entirely within the coils

maps = zoidberg.make_maps(rectangle, magnetic_field)

# Pick a poloidal slice and the next slice
yslice = 0
pol, ycoord = rectangle.getPoloidalGrid(yslice)
pol_next, ycoord_next = rectangle.getPoloidalGrid(yslice+1)

# Plot the grid points at this poloidal slice
plt.plot(pol.R, pol.Z, 'x')

# Get the coordinates which the forward map corresponds to
R_next, Z_next = pol_next.getCoordinate( maps['forward_xt_prime'][:,yslice,:], maps['forward_zt_prime'][:,yslice,:] )

plt.plot(R_next, Z_next, 'o')

plt.show()

    
