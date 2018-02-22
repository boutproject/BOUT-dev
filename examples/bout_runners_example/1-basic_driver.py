#!/usr/bin/env python

"""Driver which runs 3d_diffusion with the options given in BOUT.inp"""

from bout_runners import basic_runner

# Create the instance
my_runs = basic_runner()

# Do the run
my_runs.execute_runs()
