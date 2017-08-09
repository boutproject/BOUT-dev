#!/usr/bin/env python

"""Driver which runs 3D_diffusion using grid files."""

from bout_runners import basic_runner
from pre_and_post_processing.grid_generator import generate_grid
import os

# Generate a grid
file_name = os.path.join("grid_files","3D_diffusion_grid.nc")
generate_grid(file_name = file_name,\
              inp_path  = "data")

my_runs = basic_runner(\
            grid_file = file_name,\
            # Copy the grid file
            cpy_grid = True,\
            # Set the flag in 3D_diffusion that a grid file will be
            # used
            additional = ('flags', 'use_grid', 'true')\
            )

my_runs.execute_runs()
