
import numpy as np

from . import rzline
from . import poloidal_grid

def test_out_of_domain():
    # Create inner and outer boundaries
    inner = rzline.circle(R0=0.0, r=1.0, n=20)
    outer = rzline.circle(R0=0.0, r=2.0, n=20)

    grid = poloidal_grid.grid_elliptic(inner, outer, 10, 10)

    # Test a point in domain

    x,y = grid.findIndex(1.0, 1.0)
    r,z = grid.getCoordinate(x,y)
    assert np.allclose(r, 1.0) and np.allclose(z, 1.0)  # Should get original R,Z

    # Test out of domain

    #x,y = grid.findIndex(2.5,0)
    #print(x,y)
