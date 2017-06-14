#!/usr/bin/env python

"""Driver which runs 3D_diffusion by submitting a job to a PBS using
additional options."""

from bout_runners import PBS_runner

my_runs = PBS_runner(\
            # Although nproc is a member of basic_runner, it is used
            # together with BOUT_nodes and BOUT_ppn
            nproc                 = 4,\
            # Number of nodes to be used on the cluster
            BOUT_nodes            = 1,\
            # Specifying processor per node
            BOUT_ppn              = 4,\
            # The maximum walltime of the run
            BOUT_walltime         = '0:15:00',\
            # Specify the queue to submit to (if any)
            BOUT_queue            = None,\
            # Specify a mail to be noticed when the run has finished
            BOUT_mail             = None\
            )

# Put this in the post-processing function
my_runs.execute_runs(remove_old = True)
