#!/usr/bin/env python

from bout_runners.bout_runners import basic_runner

"""Driver which runs 3d_diffusion with other options than give in BOUT.inp"""

my_runs = basic_runner(\
            ##nproc      = 1,\
            ##directory  = 'data',\
            #solver     = ['cvode', 'ida'],\
            mms        = None,\
            atol       = None,\
            rtol       = None,\
            mxstep     = None,\
            grid_file  = None,\
            #nx         = (4, 8),\
            #ny         = [6, 7],\
            #nz         = (7, 5),\
            zperiod    = None,\
            zmin       = None,\
            zmax       = None,\
            dx         = None,\
            dy         = None,\
            dz         = None,\
            #MXG        = 1,\
            ##MYG        = 1,\
            ixseps1    = None,\
            ixseps2    = None,\
            jyseps1_1  = None,\
            jyseps1_2  = None,\
            jyseps2_1  = None,\
            jyseps2_2  = None,\
            symGlobX   = None,\
            symGlobY   = None,\
            ##ddx_first  = 'C2',\
            ##ddx_second = 'C2',\
            ##ddx_upwind = 'U1',\
            ##ddx_flux   = 'SPLIT',\
            ##ddy_first  = 'C2',\
            ##ddy_second = 'C2',\
            ##ddy_upwind = 'U1',\
            ##ddy_flux   = 'SPLIT',\
            ##ddz_first  = 'FFT',\
            #ddz_second = ['FFT', 'C2'],\
            ##ddz_upwind = 'U4',\
            ##ddz_flux   = 'SPLIT',\
            #nout       = 1,\
            #timestep   = 0.1,\
            ## additional = ('test1','ok1',1),\
            #additional = [('test1','ok1',(1, 2)), ('test2','ok2',1)],\
            sort_by    = None,\
            make       = None,\
            #restart    = 'overwrite',\
            #cpy_source = True,\
                    )

my_runs.execute_runs(remove_old = True)
