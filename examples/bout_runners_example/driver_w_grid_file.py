from bout_runners.bout_runners import basic_runner

my_instance = basic_runner(\
                    nproc      = 1,\
                    directory  = 'data',\
                    solver     = ['cvode', 'ida'],\
                    grid_file  = ['conduct_grid.nc', 'lol.nc'],\
                    ddx_first  = ['C2', 'C4'],\
                    ddx_second = 'W3',\
                    ddx_upwind = ['U1', 'W3'],\
                    ddx_flux   = ['SPLIT', 'NND'],\
                    ddy_first  = 'C2',\
                    ddy_second = ['C4', 'W3'],\
                    ddy_upwind = 'W3',\
                    ddy_flux   = 'NND',\
                    ddz_first  = 'FFT',\
                    ddz_second = ['FFT', 'C2'],\
                    #ddz_upwind = 'U4',\
                    ddz_flux   = 'SPLIT',\
                    nout       = [0, 1],\
                    timestep   = [1, 2],\
                    MXG        = 1,\
                    MYG        = 1,\
                    # additional = ('test1','ok1',1),\
                    additional = [('test1','ok1',[1, 2]), ('test2','ok2',1)],\
                    restart    = 'overwrite',\
                    cpy_source = True)

my_instance.run()
