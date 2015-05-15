from __future__ import print_function
##################################################
#            Data utilities package
#
# Generic routines, useful for all data
##################################################

import sys

print("Loading data utilities")

# Load routines from separate files
#try:
#    from plotdata import plotdata
#except:
#    print "No plotdata"
try:
    from boututils.datafile import DataFile
except:
    print("No datafile")

try:
    from boututils.file_import import file_import
except:
    print("No file_import")

try:
    from boututils.calculus import deriv, integrate
except:
    print("No calculus")

try:
    from boututils.linear_regression import linear_regression
except:
    print("No linear regression")

try:
    from boututils.shell import shell
except:
    print("No shell commands")

try:
    from boututils.ncpus import determineNumberOfCPUs
except:
    print("No determineNumberOfCPUs")

try:
    from boututils.launch import launch
except:
    print("No launch command")

try:
    from boututils.getmpirun import getmpirun
except:
    print("No getmpirun command")

try:
    from boututils.fft_integrate import fft_integrate
except:
    print("No fft_integrate command")

try:
    from boututils.mode_structure import mode_structure
except:
    print("No mode_structure command")

try:
    if sys.version_info[0]==3:
        print("polplotslice uses the VTK library through mayavi, which"+\
              " is currently only available in python 2")
    else:
        from boututils.plotpolslice import plotpolslice
except:
    print("No plotpolslice command")

try:
    from boututils.moment_xyzt import moment_xyzt
except:
    print("No moment_xyzt command")

try:
    from boututils.volume_integral import volume_integral
except:
    print("No volume_integral command")

try:
    from boututils.surface_average import surface_average
except:
    print("No surface_average command")

try:
    from boututils.showdata import showdata
except:
    print("No showdata")

try:
    from boututils.shiftz import shiftZ
except:
    print("No shiftZ")

try:
    from boututils.closest_line import closest_line
except:
    print("No closest_line")

try:
    from boututils.fft_deriv import fft_deriv
except:
    print("No fft_deriv")

try:
    from boututils.int_func  import int_func
except:
    print("No int_func")

try:
    from boututils.surface_average import surface_average
except:
    print("No surface_average ")

try:
    from boututils.efit_analyzer import View2D
except:
    print("No View2D ")

try:
    if sys.version_info[0]==3:
        print("mlab uses the VTK library through mayavi, which"+\
              " is currently only available in python 2")
    else:
        from mayavi import mlab
except:
    print("No mlab")

try:
    if sys.version_info[0]==3:
        print("anim uses the VTK library through mayavi, which"+\
              " is currently only available in python 2")
    else:
        from boututils.anim import anim
except:
    print("No anim")

try:
    if sys.version_info[0]==3:
        print("View3D uses the VTK library through mayavi, which"+\
              " is currently only available in python 2")
    else:
        from boututils.View3D import View3D
except:
    print("No View3D")
