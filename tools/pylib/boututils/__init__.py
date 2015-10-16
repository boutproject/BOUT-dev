""" Generic routines, useful for all data """

import sys

# Modules to be imported independent of version
for_all_versions = [\
                    'anim',\
                    'calculus',\
                    'closest_line',\
                    'datafile',\
                    'efit_analyzer',\
                    'fft_deriv',\
                    'fft_integrate',\
                    'file_import',\
                    'getmpirun',\
                    'int_func',\
                    'launch',\
                    'linear_regression',\
                    'mode_structure',\
                    'moment_xyzt',\
                    'ncpus',\
                    'shell',\
                    'showdata',\
                    'surface_average',\
                    'volume_integral',\
                    ]

# Check the current python version
if sys.version_info[0]>=3:
    do_import = for_all_versions
    __all__ = do_import
else:
    do_import = for_all_versions
    do_import.append('plotpolslice')
    do_import.append('View3D')
    __all__ = do_import
