#!/usr/bin/env python

"""Driver which resizes the grid after restart"""

from pre_and_post_processing.post_processing_show_the_data import show_the_data
from bout_runners import basic_runner

# Initial run
# ===========================================================================
init_run = basic_runner(nz = 8)

dmp_folder, _ =\
        init_run.execute_runs(\
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
# ===========================================================================


# Restart the run after resizing the grid
# ===========================================================================
restart_run = basic_runner(restart      = "overwrite"  ,\
                           restart_from = dmp_folder[0],\
                           nx           = 22           ,\
                           ny           = 22           ,\
                           nz           = 16           ,\
                           )

restart_run.execute_runs(\
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
# ===========================================================================
