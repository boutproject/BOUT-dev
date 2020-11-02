""" Generic routines, useful for all data """

import sys

try:
    from builtins import str
except ImportError:
    raise ImportError("Please install the future module to use Python 2")

# Modules to be imported independent of version
for_all_versions = [\
                    'calculus',\
                    'closest_line',\
                    'datafile',\
                    # 'efit_analyzer',\ # bunch pkg required
                    'fft_deriv',\
                    'fft_integrate',\
                    'file_import',\
                    'int_func',\
                    'linear_regression',\
                    'mode_structure',\
                    # 'moment_xyzt',\   # bunch pkg requried
                    'run_wrapper',\
                    'shell',\
                    'showdata',\
                    # 'surface_average',\
                    # 'volume_integral',\ #bunch pkg required
                    ]

# Check the current python version
if sys.version_info[0]>=3:
    do_import = for_all_versions
    __all__ = do_import
else:
    do_import = for_all_versions
    do_import.append('anim')
    do_import.append('plotpolslice')
    do_import.append('View3D')
    __all__ = do_import

__version__ = '0.1.4'
__name__ = 'boututils'
