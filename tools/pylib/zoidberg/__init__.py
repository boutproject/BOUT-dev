try:
    from builtins import str
except ImportError:
    raise ImportError("Please install the future module to use Python 2")

import zoidberg
from . import grid
from . import field
from . import fieldtracer
from . import plot

from .zoidberg import make_maps, write_maps

__all__ = [
    "field",
    "fieldtracer",
    "grid",
    "make_maps",
    "plot",
    "write_maps",
    "zoidberg",
]
