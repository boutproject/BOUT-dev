
import numpy as np
from . import zoidberg, grid, field

def test_make_maps_slab():
    nx = 5
    ny = 6
    nz = 7

    # Create a straight magnetic field in a slab
    straight_field = field.Slab(By=1.0, Bz=0.0, Bzprime=0.0)
    
    # Create a rectangular grid in (x,y,z)
    rectangle = grid.rectangular_grid(nx,ny,nz)

    # Calculate forwards and backwards maps
    maps = zoidberg.make_maps(rectangle, straight_field)
    
    # Check that maps has the required forward and backward index variables
    for var in ['forward_xt_prime', 'forward_zt_prime', 'backward_xt_prime', 'backward_zt_prime']:
        assert var in maps

    # Each map should have the same shape as the grid
    assert maps['forward_xt_prime'].shape == (nx,ny,nz)
    assert maps['backward_xt_prime'].shape == (nx,ny,nz)
    assert maps['forward_zt_prime'].shape == (nx,ny,nz)
    assert maps['backward_zt_prime'].shape == (nx,ny,nz)
    
    # Since this is a straight magnetic field in a simple rectangle, 
    # all the maps should be the same, and should be the identity

    identity_map_x, identity_map_z = np.meshgrid(np.arange(nx), np.arange(nz), indexing='ij')

    for y in range(ny-1):
        assert np.allclose(maps['forward_xt_prime'][:,y,:], identity_map_x)
        assert np.allclose(maps['forward_zt_prime'][:,y,:], identity_map_z)

    for y in range(1,ny):
        assert np.allclose(maps['backward_xt_prime'][:,y,:], identity_map_x)
        assert np.allclose(maps['backward_zt_prime'][:,y,:], identity_map_z)

    # The last forward map should hit a boundary
    assert np.all(maps['forward_xt_prime'][:,-1,:] < 0.0)
    assert np.all(maps['forward_zt_prime'][:,-1,:] < 0.0)

    # First backward map hits boundary
    assert np.all(maps['backward_xt_prime'][:,0,:] < 0.0)
    assert np.all(maps['backward_zt_prime'][:,0,:] < 0.0)


def test_make_maps_straight_stellarator():
    nx = 5
    ny = 6
    nz = 7
    
    # Create magnetic field
    magnetic_field = field.StraightStellarator(radius = np.sqrt(2.0))
    
    # Create a rectangular grid in (x,y,z)
    rectangle = grid.rectangular_grid(nx,ny,nz,
                                      Lx = 1.0, Lz = 1.0, Ly = 10.0,
                                      yperiodic = True)
    
    # Here both the field and and grid are centred at (x,z) = (0,0)
    # and the rectangular grid here fits entirely within the coils

    maps = zoidberg.make_maps(rectangle, magnetic_field)

    
