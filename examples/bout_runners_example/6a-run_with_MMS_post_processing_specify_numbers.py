#!/usr/bin/env python

"""Driver which runs 3D_diffusion and performs a MMS test by specifying
the grids by hand (see 6b-run_with_MMS_post_processing_grid_file.py to
see how to do the same is done by using grid files)"""

from bout_runners import basic_runner
from pre_and_post_processing.post_processing_MMS import perform_MMS_test

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
            nx = (5, 8, 16),\
            ny = (5, 8, 16),\
            nz = (4, 8, 16),\
            # Additional (put here to illustrate the sorting)
            series_add = (('cst','D_par',(1,2)), ('cst','D_perp',(0.5,1))),\
            # Since we would like to do a MMS test, we would like to run
            # the runs in a particular order. In this example, we would
            # like to run all the possible spatial variables before
            # doing the test. Hence we would like the spatial domain
            # option to be the fastest varying.
            # Since we have put post_process_after_every_run = False in
            # the run function below, the processing function being
            # called when all possibilities of the fastest variable has
            # been run.
            sort_by = 'spatial_domain'\
            # Some additional sorting examples:
            #
            # This returns an error, stating the sorting possibilities
            # (which will depend on the member data of this object)
            # sort_by = 'uncomment_me'\
            #
            # In this example cst:D_par will be the fastest varying
            # variable, followed by the spatial_domain. The post
            # processing function will be called when all possibilities
            # of these variables has been run
            # sort_by = ('cst:D_par', 'spatial_domain')\
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
