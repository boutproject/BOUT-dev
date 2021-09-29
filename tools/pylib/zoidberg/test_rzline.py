from . import rzline

import numpy as np


def test_circular_boundaries():
    inner = rzline.circle(R0=1.0, r=1.0, n=20)
    outer = rzline.circle(R0=1.0, r=2.0, n=20)

    assert np.amax(inner.R) <= 2.0
    assert np.amin(inner.R) >= 0.0


def test_distance():
    # Check the RZline.distance() function

    n = 50
    inner = rzline.circle(R0=1.0, r=1.0, n=n)
    outer = rzline.circle(R0=1.0, r=2.0, n=n)

    dist = inner.distance()

    # Check size of the array
    assert len(dist) == n + 1

    # Check start and end values
    assert np.allclose(dist[0], 0.0)
    assert np.allclose(dist[-1], 2.0 * np.pi)

    # Check that dist is monotonically increasing
    assert np.all(np.greater(dist[1:] - dist[:-1], 0.0))

    dist_outer = outer.distance()
    assert np.allclose(dist_outer[-1], 4.0 * np.pi)


def test_order_by_distance():
    # Check the RZline.equallySpaced function
    inner = rzline.circle(R0=1.0, r=1.0, n=20)
    outer = rzline.circle(R0=1.0, r=2.0, n=20)

    new_inner = inner.equallySpaced(n=10)

    assert len(new_inner.R) == 10


def test_line_from_points():
    # Create a shaped periodic line
    original = rzline.shaped_line(
        R0=3.0, a=1.0, elong=1.0, triang=0.4, indent=1.0, n=20
    )

    # Permute points
    from numpy import random

    random.seed(1235)  # Fix seed so same every time
    indx = random.permutation(original.R.size)

    R = original.R[indx]
    Z = original.Z[indx]

    reconstructed = rzline.line_from_points(R, Z)

    # Reconstructed should now be a rotated version of original
    match = False
    for i in range(len(R)):
        if np.allclose(original.R, np.roll(reconstructed.R, i)) and np.allclose(
            original.Z, np.roll(reconstructed.Z, i)
        ):
            match = True
            break

    assert match
