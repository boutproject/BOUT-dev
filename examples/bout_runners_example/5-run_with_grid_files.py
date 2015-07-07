#!/usr/bin/env python

"""Driver which runs 3D_diffusion using grid files."""

from bout_runners.bout_runners import basic_runner
from grid_generator import generate_grid

# Generate a grid
generate_grid(file_name = "3D_diffusion_grid",\
              inp_path  = "data")

my_runs = basic_runner(\
            grid_file = "3D_diffusion_grid.nc",\
            # Copy the grid file
            cpy_grid = True,\
            # Set the flag in 3D_diffusion that a grid file will be
            # used
            additional = ('flags', 'use_grid', 'true')\
            )

my_runs.execute_runs()
