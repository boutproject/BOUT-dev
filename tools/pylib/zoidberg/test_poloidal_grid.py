from . import poloidal_grid as grid

import numpy as np

def test_circular_boundaries():
    result = grid.circular_boundaries(R0=1.0, rin=1.0, rout=2.0, n=20)
    
    assert len(result) == 2 # Should return two values
    
    inner, outer = result
    
    assert np.amax(inner.R) <= 2.0
    assert np.amin(inner.R) >= 0.0
    
    
def test_distance():
    # Check the RZline.distance() function
    
    n = 50
    inner,outer = grid.circular_boundaries(R0=1.0, rin=1.0, rout=2.0, n=n)
    
    dist = inner.distance()

    # Check size of the array
    assert len(dist) == n+1
    
    # Check start and end values
    assert np.allclose(dist[0], 0.0)
    assert np.allclose(dist[-1], 2.*np.pi)
    
    # Check that dist is monotonically increasing
    assert np.all(np.greater(dist[1:] - dist[:-1], 0.0)) 
    
    dist_outer = outer.distance()
    assert np.allclose(dist_outer[-1], 4.*np.pi)
    
def test_order_by_distance():
    # Check the RZline.orderByDistance function
    
    inner,outer = grid.circular_boundaries(R0=1.0, rin=1.0, rout=2.0, n=20)
    
    new_inner = inner.orderByDistance(n=10)
    
    assert len(new_inner.R) == 10
    
    
def test_grid_annulus():
    
    inner,outer = grid.circular_boundaries(R0=1.0, rin=1.0, rout=2.0, n=50)
    
    grid.grid_annulus(inner, outer, 5, 10, show=True)
    
