#!/usr/bin/env python

"""Driver which restarts a scan, given a restart function"""

from bout_runners import PBS_runner
from pre_and_post_processing.restart_from_func import restart_from_func

scan = (("cst", "D_perp", (1.0,5.5)),\
        ("cst", "D_par",  (1.5,2.5))
       )

# Initial runs
# ===========================================================================
init_run = PBS_runner(additional = scan)

dmp_folder, PBS_ids =\
        init_run.execute_runs()

one_of_the_restart_paths_in_scan = dmp_folder[0]
# ===========================================================================


# Restart the scan
# ===========================================================================
restart_run = PBS_runner(nout         = 5                ,\
                         restart      = "overwrite"      ,\
                         restart_from = restart_from_func,\
                         additional   = scan             ,\
                        )

restart_run.execute_runs(\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
                         # Below are the kwargs given to the
                         # restart_from_func
                         one_of_the_restart_paths_in_scan =\
                         one_of_the_restart_paths_in_scan,\
                         scan_parameters = ("D_perp", "D_par"),\
                        )
# ===========================================================================
