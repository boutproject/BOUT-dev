import zoidberg
from . import grid
from . import field
from . import fieldtracer
try:
    from . import plot
except:
    #no plotting available
    pass

from .zoidberg import make_maps, write_maps
