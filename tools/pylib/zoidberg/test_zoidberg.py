
from . import zoidberg, grid, field

def test_make_maps_slab():

    # Create a straight magnetic field in a slab
    straight_field = field.Slab(By=1.0, Bz=0.0, Bzprime=0.0)
    
    # Create a rectangular grid in (x,y,z)
    rectangle = grid.rectangular_grid(10,10,10)

    # Calculate forwards and backwards maps
    maps = zoidberg.make_maps(rectangle, straight_field)

    print maps
