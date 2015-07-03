#!/usr/bin/env python

"""Driver which runs 3d_diffusion and performs a MMS test by specifying
the grids by hand (see 6b-run_with_MMS_post_processing_grid_file.py to
see how to do the same using grid files)"""

from bout_runners.bout_runners import basic_runner
from post_processing_MMS import perform_MMS_test

# Specifying the grids
# For xl
# nx = [8,16,32]
# ny = [4,4,4]
# nz = [1,1,1]
# For yl
nx = [4,4,4]
ny = [8,16,32]
nz = [1,1,1]

my_runs = basic_runner(\
            # TODO: del
            solver = 'rk4',\
            nproc = 4,\
            # Set the directory
            directory = 'MMS',\
            # Set the time domain
            nout     = 1,\
            timestep = 1,\
            # Set mms to true
            mms = True,\
            # Set the spatial domain
            nx = nx,\
            ny = ny,\
            nz = nz,\
            # TODO: del
            make = True,\
            # TODO: del
            # additional = ('solver', 'monitor_timestep', 'true'),\
            )

# Put this in the post-processing function
my_runs.execute_runs(\
                     # TODO: del
                     remove_old = True,\
                     post_processing_function = perform_MMS_test,\
                     # As we need several runs in order to perform the
                     # MMS test, this needs to be false
                     post_process_after_every_run = False,\
                     # Below are the kwargs arguments being passed to
                     # perform_MMS_test
                     extension = 'png',\
                     show_plot = True
                    )
