##################################################
#            BOUT++ data package
#
# Routines for exchanging data to/from BOUT++
#
##################################################

print "Loading BOUT++ data routines"

# Load routines from separate files
try:
    from collect import collect
except:
    print "Sorry, no collect"

try:
    from pol_slice import pol_slice
except:
    print "Sorry, no pol_slice"


