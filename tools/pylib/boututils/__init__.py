##################################################
#            Data utilities package
#
# Generic routines, useful for all data
##################################################

print "Loading data utilities"

# Load routines from separate files
try:
    from showdata import showdata
except:
    print "No showdata"

try:
    from plotdata import plotdata
except:
    print "No plotdata"

try:
    from datafile import DataFile
except:
    print "No datafile"

try:
    from calculus import deriv, integrate
except:
    print "No calculus"
