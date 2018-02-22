#!/usr/bin/env python

import numpy as np
import zoidberg
   
field = zoidberg.field.GEQDSK("g014220.00200") # Read magnetic field

grid = zoidberg.grid.rectangular_grid(100, 10, 100,
                                      1.5-0.1, # Range in R (max - min)
                                      2*np.pi, # Toroidal angle
                                      3., # Range in Z
                                      xcentre=(1.5+0.1)/2, # Middle of grid in R
                                      yperiodic=True) # Periodic in toroidal angle

# Create the forward and backward maps
maps = zoidberg.make_maps(grid, field)

# Save to file
zoidberg.write_maps(grid, field, maps, gridfile="tokamak.fci.nc")

#############################################################################
# Plot maps

zoidberg.plot.plot_forward_map(grid, maps)

