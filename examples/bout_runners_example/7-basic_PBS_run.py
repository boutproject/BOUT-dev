#!/usr/bin/env python

"""
Driver which runs 3D_diffusion by submitting a job to a Portable Batch System
(PBS)
"""

from bout_runners import PBS_runner

my_runs = PBS_runner()

my_runs.execute_runs()
