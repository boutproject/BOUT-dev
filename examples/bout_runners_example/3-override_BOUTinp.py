#!/usr/bin/env python

"""Driver which runs 3d_diffusion with other options than given in BOUT.inp"""

from bout_runners import basic_runner

my_runs = basic_runner(\
            # Number of processors
            nproc      = 2,\
            # Directory of the inp file
            directory  = 'data',\
            # Set the solver option
            solver     = 'rk4',\
            mms        = False,\
            atol       = 1.0e-8,\
            rtol       = 1.0e-8,\
            mxstep     = 10000000,\
            # Spatial domain option
            nx         = 19,\
            ny         = 17,\
            nz         = 16,\
            # These can be set if needed
            zperiod    = None,\
            zmin       = None,\
            zmax       = None,\
            # These are not set here, but the code handles them
            # internally
            dx         = None,\
            dy         = None,\
            dz         = None,\
            # The same as in BOUT.inp
            # (Setting them to a different value doesn't make much sense)
            MXG        = 1,\
            MYG        = 1,\
            # These can also be set
            ixseps1    = None,\
            ixseps2    = None,\
            jyseps1_1  = None,\
            jyseps1_2  = None,\
            jyseps2_1  = None,\
            jyseps2_2  = None,\
            symGlobX   = None,\
            symGlobY   = None,\
            # The differencing option
            ddx_first  = 'C2',\
            ddx_second = 'C2',\
            ddx_upwind = 'U1',\
            ddx_flux   = 'SPLIT',\
            ddy_first  = 'C2',\
            ddy_second = 'C2',\
            ddy_upwind = 'U1',\
            ddy_flux   = 'SPLIT',\
            ddz_first  = 'FFT',\
            ddz_second = 'FFT',\
            ddz_upwind = 'U4',\
            ddz_flux   = 'SPLIT',\
            # Temporal domain option
            nout       = 11,\
            timestep   = 0.02,\
            # Additional options
            # (An example ofadditional options run in series is found in
            #  6a-run_with_MMS_post_processing_specify_numbers.py)
            # tuple[0] - section name
            # tuple[1] - variable name for the section
            # tuple[2] - value of the variable name in the section
            additional = (('cst','D_perp',5), ('cst', 'D_par', 0.5)),\
            # Can set this to overwrite or append
            restart    = None,\
            # Will copy the source file
            cpy_source = True,\
            # Will remake the file
            make       = True,\
            # Code will return an error if False, due to the mismatch
            # between nx, ny and nproc
            allow_size_modification = True)

my_runs.execute_runs(\
        # Remove eventually old data
        remove_old = True\
        )
