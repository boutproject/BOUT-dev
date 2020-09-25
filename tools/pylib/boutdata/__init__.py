""" Routines for exchanging data to/from BOUT++ """

try:
    from builtins import str
except ImportError:
    raise ImportError("Please install the future module to use Python 2")

# Import this, as this almost always used when calling this package
from boutdata.collect import collect, attributes

__all__ = ["attributes", "collect", "gen_surface", "pol_slice"]

__version__ = '0.1.2'
__name__ = 'boutdata'
