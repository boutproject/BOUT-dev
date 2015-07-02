#!/usr/bin/env python

"""Driver which runs 3d_diffusion and performs a MMS test."""

from bout_runners.bout_runners import basic_runner
from grid_generator import generate_grid
from post_processing_MMS import MMS_test

# TODO: Make a a version and b version of this file, do not let user
# specify

# The test can be done by using grid files, or by specifying nx, ny and
# nz. The choice is yours

# Options
# 'grid_file'
# 'specify_numbers'
MMS_by = 'grid_file'

# Just checking if you managed to type MMS by correctly
if MMS_by != 'grid_file' or MMS_by != 'specify_numbers':
    raise TypeError("MMS must be 'grid_file' or 'specify_numbers'")

# Do MMS test by the grid file
if MMS_by == 'grid_file':

# Do MMS test by specifying nx, ny and nz
elif MMS_by == 'specify_numbers':

# Generate a grid
generate_grid(file_name="3D_diffusion_grid")

my_runs = basic_runner(\
            grid_file = "3D_diffusion_grid.nc",\
            # Copy the grid file
            cpy_grid = True\
            )

my_runs.execute_runs()
