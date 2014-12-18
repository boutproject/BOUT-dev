##################################################
#            Data utilities package
#
# Generic routines, useful for all data
##################################################
import traceback

print("Loading data utilities")

# Load routines from separate files
#try:
#    from plotdata import plotdata
#except:
#    print "No plotdata"
try:
    from datafile import DataFile
except:
    print("No datafile")

try:
    from file_import import file_import
except:
    print("No file_import")

try:
    from calculus import integrate, deriv 
except:
    print("No calculus")

try:
    from linear_regression import linear_regression
except:
    print("No linear regression")

try:
    from shell import shell
except:
    print("No shell commands")

try:
    from ncpus import determineNumberOfCPUs
except:
    print("No determineNumberOfCPUs")

try:
    from launch import launch
except:
    print("No launch command")

try:
    from getmpirun import getmpirun
except:
    print("No getmpirun command")
 
try:
    from fft_integrate import fft_integrate
except:
    print "No fft_integrate command"

try:
    from mode_structure import mode_structure
except:
    print "No mode_structure command"

try:
    from plotpolslice import plotpolslice
except:
    print "No plotpolslice command"

try:
    from moment_xyzt import moment_xyzt
except:
    print "No moment_xyzt command"

try:
    from volume_integral import volume_integral
except:
    print "No volume_integral "

try:
    from closest_line import closest_line
except:
    print "No closest_line"

try:
    from fft_deriv import fft_deriv
except:
    print "No fft_deriv"
 
try:
    from int_func  import int_func  
except:
    print "No int_func"
 
try:
    from surface_average import surface_average
except:
    print "No surface_average "

try:
    from efit_analyzer import View2D
except:
    print "No View2D "

try:
    from mayavi import mlab
except:
    print("No mlab")

try:
    from anim import anim
except:
    print traceback.format_exc()


try:
    from View3D import View3D
except:
    print traceback.format_exc()

