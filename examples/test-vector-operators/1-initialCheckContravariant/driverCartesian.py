#!/usr/bin/env python

"""Driver which checks cartesian perturbation."""

from bout_runners.bout_runners import basic_runner
import numpy as np
from pythonRoutines.vecFieldPlotCart import plotVecField
# The options for the run
# =============================================================================
# General options
restart    = None
remove_old = True
directory  = "cartesian"
# Shall we make?
make       = True
# The number of processors
nproc = 4
# =============================================================================


# The options for the post processing function
# =============================================================================
scale = 50
# =============================================================================


# Create the runner
# =============================================================================
my_runs = basic_runner(\
            directory  = directory ,\
            nproc      = nproc ,\
            # Copy the source file
            cpy_source = True  ,\
            make       = make  ,\
            restart    = restart,\
            )
# =============================================================================


# Perform the run
# =============================================================================
my_runs.execute_runs(\
        remove_old = remove_old,\
        post_processing_function = plotVecField,\
        # This function will be called every time after
        # performing a run
        post_process_after_every_run = True,\
        # Below are the kwargs arguments being passed to
        # the post processing function
        scale = scale,\
        )
# =============================================================================
