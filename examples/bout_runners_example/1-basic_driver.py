#!/usr/bin/env python

"""Driver which runs conduction with the options given in BOUT.inp"""

from bout_runners.bout_runners import basic_runner

# Create the instance
my_run = basic_runner(make=True)

# Do the run
my_run.run(remove_old=True)
