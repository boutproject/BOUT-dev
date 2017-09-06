#!/usr/bin/env python

"""Driver which restarts a scan, given a restart function"""

from pre_and_post_processing.post_processing_show_the_data import show_the_data
from pre_and_post_processing.restart_from_func import restart_from_func
from bout_runners import basic_runner

scan = (("cst", "D_perp", (1.0,5.5)),\
        ("cst", "D_par",  (1.5,2.5))
       )

# Given that the runs has already been performed
only_post_process = False

# Initial runs
# ===========================================================================
init_run = basic_runner(additional = scan)

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

one_of_the_restart_paths_in_scan = dmp_folder[0]
# ===========================================================================


# Restart the scan
# ===========================================================================
if only_post_process:
    restart = None
else:
    restart = "overwrite"

restart_run = basic_runner(nout         = 5                ,\
                           restart      = restart          ,\
                           restart_from = restart_from_func,\
                           additional   = scan             ,\
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
                         z = slice(0,None),\
                         # Below are the kwargs given to the
                         # restart_from_func
                         one_of_the_restart_paths_in_scan =\
                         one_of_the_restart_paths_in_scan,\
                         scan_parameters = ["D_perp", "D_par"],\
                        )
# ===========================================================================
