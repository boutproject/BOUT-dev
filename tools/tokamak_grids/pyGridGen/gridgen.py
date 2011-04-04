# Grid generation code in Python
#
# Takes an RZ equilibrium from a variety of sources (e.g. EFIT),
# and produces a set of flux surfaces inside and outside separatrix.
#
# This can then be put into the grid2bout code
# 
# Ben Dudson, 2009

try:
    import numpy as np
except ImportError:
    print "ERROR: NumPy module not available"
    raise

def gridgen(psi_rz):
    



