##################################################
#            BOUT++ data package
#
# Routines for exchanging data to/from BOUT++
#
##################################################

print("Loading BOUT++ data routines")

# Load routines from separate files
from boutdata.collect import collect

from boutdata.pol_slice import pol_slice

from boutdata.gen_surface import gen_surface

