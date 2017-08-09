#!/usr/bin/env python

"""Driver which runs 3D_diffusion and performs a MMS test by specifying
the grids by using grid_files (see
6a-run_with_MMS_post_processing_specify_numbers.py to see how to do the
same is done by specifying the grid manually)"""

from bout_runners import basic_runner
from pre_and_post_processing.post_processing_MMS import perform_MMS_test
from pre_and_post_processing.grid_generator import  generate_grid
import os

# Generate the grids
# Specify the grid dimensions
grid_numbers = (5, 8, 16)
# Make an append able list
grid_files = []
for grid_number in grid_numbers:
    file_name = os.path.join("grid_files","grid_file_{}.nc".format(grid_number))
    # Generate the grids
    generate_grid(nx        = grid_number,\
                  ny        = grid_number,\
                  nz        = grid_number,\
                  inp_path  = 'MMS'      ,\
                  file_name = file_name)
    # Append the grid_files list
    grid_files.append(file_name)

my_runs = basic_runner(\
            nproc     = 1,\
            # Set the directory
            directory  = 'MMS',\
            # Set the time domain
            nout       = 1,\
            timestep   = 1,\
            # Set mms to true
            mms        = True,\
            # Set the spatial domain
            grid_file  = grid_files,\
            # Set the flag in 3D_diffusion that a grid file will be
            # used
            additional = ('flags','use_grid','true'),\
            # Copy the grid file
            cpy_grid   = True,\
            # Sort the runs by the spatial domain
            sort_by    = 'grid_file'
            )

# Put this in the post-processing function
my_runs.execute_runs(\
                     post_processing_function = perform_MMS_test,\
                     # As we need several runs in order to perform the
                     # MMS test, this needs to be false
                     post_process_after_every_run = False,\
                     # Below are the kwargs arguments being passed to
                     # perform_MMS_test
                     extension = 'png',\
                     show_plot = True\
                    )
