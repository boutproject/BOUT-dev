import numpy as np

from . import fieldtracer
from . import field


def test_slab():
    mag_field = field.Slab(By=0.5, Bz=0.1, Bzprime=0.0)

    tracer = fieldtracer.FieldTracer(mag_field)

    coords = tracer.follow_field_lines(0.2, 0.3, [1.5, 2.5])

    assert coords.shape == (2, 1, 2)
    assert np.allclose(coords[:, 0, 0], 0.2)  # X coordinate

    assert np.allclose(coords[0, 0, 1], 0.3)
    assert np.allclose(coords[1, 0, 1], 0.3 + 0.1 / 0.5)  # Z coordinate

    coords = tracer.follow_field_lines([0.2, 0.3, 0.4], 0.5, [1.0, 2.5])

    assert coords.shape == (2, 3, 2)
    assert np.allclose(coords[:, 0, 0], 0.2)
    assert np.allclose(coords[:, 1, 0], 0.3)
    assert np.allclose(coords[:, 2, 0], 0.4)
    assert np.allclose(coords[1, :, 1], 0.5 + 1.5 * 0.1 / 0.5)  # Z coord


def test_FieldTracerReversible_slab():
    mag_field = field.Slab(By=0.5, Bz=0.1, Bzprime=0.0)

    tracer = fieldtracer.FieldTracerReversible(mag_field)

    coords = tracer.follow_field_lines(0.2, 0.3, [1.5, 2.5])

    assert coords.shape == (2, 1, 2)
    assert np.allclose(coords[:, 0, 0], 0.2)  # X coordinate

    assert np.allclose(coords[0, 0, 1], 0.3)
    assert np.allclose(coords[1, 0, 1], 0.3 + 0.1 / 0.5)  # Z coordinate

    coords = tracer.follow_field_lines([0.2, 0.3, 0.4], 0.5, [1.0, 2.5])

    assert coords.shape == (2, 3, 2)
    assert np.allclose(coords[:, 0, 0], 0.2)
    assert np.allclose(coords[:, 1, 0], 0.3)
    assert np.allclose(coords[:, 2, 0], 0.4)
    assert np.allclose(coords[1, :, 1], 0.5 + 1.5 * 0.1 / 0.5)  # Z coord


def test_poincare():
    mag_field = field.Slab(By=0.5, Bz=0.1, Bzprime=0.0)

    result, y_slices = fieldtracer.trace_poincare(
        mag_field, 0.0, 0.0, 1.0, nplot=3, revs=5, nover=1
    )

    assert y_slices.size == 3
    assert result.shape == (5, y_slices.size, 1, 2)
