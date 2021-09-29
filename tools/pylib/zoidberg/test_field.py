import numpy as np
from .field import Slab, CurvedSlab


def test_slab():
    # Slab with no shear (Bzprime = 0)
    slab = Slab(By=0.5, Bz=0.2, Bzprime=0.0)

    # Create an array of coordinates
    x = np.arange(10)
    z = x

    Bx = slab.Bxfunc(x, z, 0.0)
    assert Bx.shape == x.shape
    assert np.allclose(Bx, 0.0)

    By = slab.Byfunc(x, z, 0.0)
    assert By.shape == x.shape
    assert np.allclose(By, 0.5)

    Bz = slab.Bzfunc(x, z, 0.0)
    assert Bz.shape == x.shape
    assert np.allclose(Bz, 0.2)


def test_curved_slab():
    # Slab with no shear (Bzprime = 0)
    slab = CurvedSlab(By=0.5, Bz=0.2, Bzprime=0.0)

    # Create an array of coordinates
    x = np.arange(10)
    z = x

    Bx = slab.Bxfunc(x, z, 0.0)
    assert Bx.shape == x.shape
    assert np.allclose(Bx, 0.0)

    By = slab.Byfunc(x, z, 0.0)
    assert By.shape == x.shape
    assert np.allclose(By, 0.5)

    Bz = slab.Bzfunc(x, z, 0.0)
    assert Bz.shape == x.shape
    assert np.allclose(Bz, 0.2)
