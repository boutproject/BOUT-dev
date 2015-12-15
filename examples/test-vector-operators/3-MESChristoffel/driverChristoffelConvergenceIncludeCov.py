#!/usr/bin/env python

"""
Driver which checks the Christoffel symbols. In this case the covariant
metric tensor elements are included in the BOUT.inp file which we run
from.
"""

from bout_runners.bout_runners import basic_runner
import numpy as np
from pythonRoutines.checkConvergence import convergenceTest
# The options for the run
# =============================================================================
directory = 'includeCov'
# The grid
ny = [4,8,16,32]
nz = [4,8,16,32]
nx = [n + 2 for n in ny]
# General options
restart    = None
remove_old = True
# Shall we make?
make       = True
# The number of processors
nproc = 4
# =============================================================================


# Create the runner
# =============================================================================
my_runs = basic_runner(\
            directory = directory,\
            # Size of the grid
            nx = nx,\
            ny = ny,\
            nz = nz,\
            nproc      = nproc ,\
            # Copy the source file
            cpy_source = True  ,\
            make       = make  ,\
            restart    = restart,\
            # Sort by the spatial domain
            sort_by = 'spatial_domain'\
            )
# =============================================================================


# Perform the run
# =============================================================================
my_runs.execute_runs(\
        remove_old = remove_old,\
        post_processing_function = convergenceTest,\
        # This function will be called every time after
        # performing a run
        post_process_after_every_run = False,\
        # Below are the kwargs arguments being passed to
        # the post processing function
        showPlot = False
        )
# =============================================================================
