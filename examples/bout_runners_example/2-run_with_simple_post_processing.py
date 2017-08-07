#!/usr/bin/env python

"""Driver which runs 3d_diffusion, and calls the function show_the_data when done"""

from pre_and_post_processing.post_processing_show_the_data import show_the_data
from bout_runners import basic_runner


my_runs = basic_runner()

# Put this in the post-processing function
my_runs.execute_runs(\
                     post_processing_function = show_the_data,\
                     # This function will be called every time after
                     # performing a run
                     post_process_after_every_run = True,\
                     # Below are the kwargs arguments being passed to
                     # show_the_data
                     t = slice(0,None),\
                     x = 1,\
                     y = slice(0,None),\
                     z = slice(0,None)\
                    )
