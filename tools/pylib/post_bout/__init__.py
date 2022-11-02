from __future__ import print_function
from __future__ import absolute_import

##################################################
#            BOUT++ data package
#
# Routines for examining simulation results for BOUT++
#
##################################################

print("Loading BOUT++ post processing routines")

# Load routines from separate files
import sys
import os

try:

    boutpath = os.environ["BOUT_TOP"]
    pylibpath = boutpath + "/tools/pylib"
    boutdatapath = pylibpath + "/boutdata"
    boututilpath = pylibpath + "/boututils"
    allpath = [boutpath, pylibpath, boutdatapath, boututilpath]
    [sys.path.append(elem) for elem in allpath]
    print(sys.path)

    # sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')

    print("in post_bout/__init__.py")

    import matplotlib

    matplotlib.use("pdf")  # savemovie must be called as a diff. sesssion

    import gobject
    import numpy as np
except ImportError:
    print("can't find the modules I need, you fail")
    sys.exit()  # no point in going on


# import some bout specific modules
try:
    import boutdata
    import boututils
except:
    print("can't find bout related modules, you fail")

# import some home-brewed modules


# create some aliases


try:
    from read_grid import read_grid
except:
    print("Sorry, no read_grid")

try:
    from .read_inp import parse_inp, read_inp, read_log, metadata
except:
    print("Sorry no parse_inp")

try:
    from .read_cxx import read_cxx, get_evolved_cxx, no_comment_cxx
except:
    print("Sorry no read_cxx")

try:
    from post_bout import save, read
except:
    print("Sorry, no show")

try:
    from .basic_info import basic_info, fft_info
except:
    print("Sorry, no basic_info")

try:
    from .pb_corral import corral, LinRes, subset
except:
    print("No corral")

try:
    from . import ListDict
except:
    print("No ListDict")

try:
    # from rotate_mp import rotate
    from rotate2 import rotate
except:
    print("No rotate")
