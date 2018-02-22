#!/usr/bin/env python

"""Driver which runs 3d_diffusion for several different combinations of
the input. """

from bout_runners import basic_runner

# With a few exceptions: All variables in the constructor can be given
# as an iterable.
# When execute_runs is called, bout_runners will run all combinations of
# the member data
my_runs = basic_runner(\
            # nx, ny and nz must be of the same size as they constitute
            # one "part" of the combination (i.e. there will be no
            # internal combination between the elements in nx, ny and
            # nz)
            nx         = (9, 18),\
            ny         = (6, 12),\
            nz         = (8, 16),\
            # nout and timestep must be of the same dimension for the
            # same reason as mention above
            nout       = (10,   11,    12),\
            timestep   = (0.01, 0.01, 0.01),\
            # The differencing option
            ddz_second = ('FFT','C2'),\
            # Additional options
            additional = (('cst','D_perp',(1, 2)))\
            )

# Execute all the runs
# 2 runs for each combination of nx, ny, nz
# 3 runs for each combination of nout and timestep
# 2 runs for each combination in ddz_second
# 2 runs for each combination of cst:const:value
# In total: 24 runs
my_runs.execute_runs()

# NOTE: If you feel that the explanation of the combinations was bad,
#       have a look at the last lines of data/run_log.txt to see what
#       runs have been performed after this run
