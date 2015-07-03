#!/usr/bin/env python

"""Driver which runs 3d_diffusion and performs a MMS test by specifying
the grids by hand (see 6b-run_with_MMS_post_processing_grid_file.py to
see how to do the same using grid files)"""

from bout_runners.bout_runners import basic_runner
from post_processing_MMS import perform_MMS_test

my_runs = basic_runner(\
            nproc = 1,\
            # Set the directory
            directory = 'MMS',\
            # Set the time domain
            nout     = 1,\
            timestep = 1,\
            # Set mms to true
            mms = True,\
            # Set the spatial domain
            nx = [4, 8, 16],\
            ny = [4, 8, 16],\
            nz = [4, 8, 16],\
            # Additional (put here to illustrate the sorting)
            series_add = [('cst','D_par',[1,2]), ('cst','D_perp',[0.5,1])]
            # TODO: The sorting
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
                     show_plot = True
                    )
