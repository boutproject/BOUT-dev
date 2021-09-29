from itertools import chain, product
import numpy as np
from . import zoidberg, grid, field


def test_make_maps_slab():
    nx = 5
    ny = 8
    nz = 7

    # Create a straight magnetic field in a slab
    straight_field = field.Slab(By=1.0, Bz=0.0, Bzprime=0.0)

    # Create a rectangular grid in (x,y,z)
    rectangle = grid.rectangular_grid(nx, ny, nz)

    # Two parallel slices in each direction
    nslice = 2

    # Calculate forwards and backwards maps
    maps = zoidberg.make_maps(rectangle, straight_field, nslice=nslice)

    # Since this is a straight magnetic field in a simple rectangle,
    # all the maps should be the same, and should be the identity
    identity_map_x, identity_map_z = np.meshgrid(
        np.arange(nx), np.arange(nz), indexing="ij"
    )

    # Check that maps has the required forward and backward index variables
    offsets = chain(range(1, nslice + 1), range(-1, -(nslice + 1), -1))
    field_line_maps = ["xt_prime", "zt_prime"]

    for field_line_map, offset in product(field_line_maps, offsets):
        var = zoidberg.parallel_slice_field_name(field_line_map, offset)
        print("Current field: ", var)
        assert var in maps

        # Each map should have the same shape as the grid
        assert maps[var].shape == (nx, ny, nz)

        # The first/last abs(offset) points are not valid, so ignore those
        interior_range = (
            range(ny - abs(offset)) if offset > 0 else range(abs(offset), ny)
        )
        # Those invalid points should be set to -1
        end_slice = slice(-1, -(offset + 1), -1) if offset > 0 else slice(0, -offset)
        identity_map = identity_map_x if "x" in var else identity_map_z

        for y in interior_range:
            assert np.allclose(maps[var][:, y, :], identity_map)

        # The end slice should hit a boundary
        assert np.allclose(maps[var][:, end_slice, :], -1.0)


def test_make_maps_straight_stellarator():
    nx = 5
    ny = 6
    nz = 7

    # Create magnetic field
    magnetic_field = field.StraightStellarator(radius=np.sqrt(2.0))

    # Create a rectangular grid in (x,y,z)
    rectangle = grid.rectangular_grid(
        nx, ny, nz, Lx=1.0, Lz=1.0, Ly=10.0, yperiodic=True
    )

    # Here both the field and and grid are centred at (x,z) = (0,0)
    # and the rectangular grid here fits entirely within the coils

    maps = zoidberg.make_maps(rectangle, magnetic_field)
