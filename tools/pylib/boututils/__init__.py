##################################################
#            Data utilities package
#
# Generic routines, useful for all data
##################################################

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
    from calculus import deriv, integrate
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
    print "No volume_integral command"

try:
    from surface_average import surface_average
except:
    print "No surface_average command"

try:
    from showdata import showdata
except:
    print "No showdata"
