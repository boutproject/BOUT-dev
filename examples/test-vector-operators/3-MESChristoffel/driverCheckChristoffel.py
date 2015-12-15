#!/usr/bin/env python

"""Driver which plots the individual Christoffel symbols."""

from bout_runners.bout_runners import basic_runner
import numpy as np
from pythonRoutines.checkChristoffel import plotChristoffel
# The options for the run
# =============================================================================
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
        post_processing_function = plotChristoffel,\
        # This function will be called every time after
        # performing a run
        post_process_after_every_run = True,\
        # Below are the kwargs arguments being passed to
        # the post processing function
        showRelErr     = False,\
        showBOUT       = True,\
        showAnalytical = False,\
        xguards        = True,\
        showOnlyThese  = ['G1_11', 'G2_13']
        )
# =============================================================================
