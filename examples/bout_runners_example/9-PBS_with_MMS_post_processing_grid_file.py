#!/usr/bin/env python

"""Driver which runs 3D_diffusion by submitting a job to a PBS and
performs a MMS test by specifying the grids by using grid_files."""

from bout_runners import PBS_runner
from pre_and_post_processing.post_processing_MMS import perform_MMS_test
from pre_and_post_processing.grid_generator import  generate_grid
import os

# Generate the grids
# Specify the grid dimensions
grid_numbers = (8, 16, 32)
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
    grid_files.append(file_name)

my_runs = PBS_runner(\
            # Specify the numbers used for the BOUT runs
            nproc                 = 4,\
            BOUT_nodes            = 1,\
            BOUT_ppn              = 4,\
            BOUT_walltime         = '0:15:00',\
            BOUT_queue            = None,\
            BOUT_mail             = None,\
            # Specify the numbers used for the post processing
            post_process_nproc    = 1,\
            post_process_nodes    = 1,\
            post_process_ppn      = 1,\
            post_process_walltime = '0:05:00',\
            post_process_queue    = None,\
            post_process_mail     = None,\
            # Set the directory
            directory             = 'MMS',\
            # Set the time domain
            nout                  = 1,\
            timestep              = 1,\
            # Set mms to true
            mms                   = True,\
            # Set the spatial domain
            grid_file             = grid_files,\
            # Set the flag in 3D_diffusion that a grid file will be
            # used
            additional            = ('flags','use_grid','true'),\
            # Add some additional option
            series_add            = (('cst','D_par' ,(1,2)),\
                                     ('cst','D_perp',(0.5,1))),\
            # Copy the grid file
            cpy_grid              = True,\
            # Sort the runs by the spatial domain
            sort_by               = 'grid_file'
            )

# Put this in the post-processing function
my_runs.execute_runs(\
                     remove_old = True,\
                     post_processing_function = perform_MMS_test,\
                     # As we need several runs in order to perform the
                     # MMS test, this needs to be false
                     post_process_after_every_run = False,\
                     # Below are the kwargs arguments being passed to
                     # perform_MMS_test
                     extension = 'png',\
                     show_plot = False\
                    )
