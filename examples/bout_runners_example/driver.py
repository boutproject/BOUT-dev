from bout_runners.bout_runners import basic_runner

my_instance = basic_runner(\
                    #nproc      = 1,\
                    #directory  = 'data',\
                    solver     = ['cvode', 'ida'],\
                    nx         = (4, 8),\
                    ny         = [6, 7],\
                    nz         = (7, 5),\
                    #ddx_first  = 'C2',\
                    #ddx_second = 'C2',\
                    #ddx_upwind = 'U1',\
                    #ddx_flux   = 'SPLIT',\
                    #ddy_first  = 'C2',\
                    #ddy_second = 'C2',\
                    #ddy_upwind = 'U1',\
                    #ddy_flux   = 'SPLIT',\
                    #ddz_first  = 'FFT',\
                    ddz_second = ['FFT', 'C2'],\
                    #ddz_upwind = 'U4',\
                    #ddz_flux   = 'SPLIT',\
                    nout       = 1,\
                    timestep   = 0.1,\
                    MXG        = 1,\
                    #MYG        = 1,\
                    # additional = ('test1','ok1',1),\
                    additional = [('test1','ok1',(1, 2)), ('test2','ok2',1)],\
                    restart    = 'overwrite',\
                    cpy_source = True)

my_instance.run()
