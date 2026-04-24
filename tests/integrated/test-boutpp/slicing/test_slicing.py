#!/usr/bin/env python3
# requires boutpp
# requires not make
import numpy as np
import numpy.testing as npt
import boutpp as bc
import inspect

def test_slicing():

    bc.init("-d test")

    mesh = bc.Mesh.getGlobal()
    field = bc.Field3D.fromMesh(mesh)
    ndat = np.random.random(field.shape)
    field[:] = ndat

    npt.assert_allclose(ndat, field[:])

    examples = [
        lambda x: x[2],
        lambda x: x[-2],
        lambda x: x[1, 3],
        lambda x: x[1, -1],
        lambda x: x[0],
        lambda x: x[0][2],
        lambda x: x[1:7:2],
        lambda x: x[-2:10],
        lambda x: x[-3:3:-1],
        lambda x: x[5:],
        lambda x: x[1:2],
        lambda x: x[..., 0],
        lambda x: x[:, :, 0],
        lambda x: x[np.array([3, 3, 1, 8])],
        lambda x: x[np.array([3, 3, -3, 8])],
        lambda x: x[np.array([1, -1])],
        lambda x: x[np.array([3, 4])],
        lambda x: x[[0, 1, 2], [0, 1, 0]],
        lambda x: x[1:2, 1:3],
        lambda x: x[1:2, [1, 2]],
        lambda x: x[ndat < 0],
    ]

    print(field.shape)
    for ex in examples:
        print("testing", inspect.getsource(ex))
        nout = ex(ndat)
        fout = ex(field)
        npt.assert_allclose(
            fout, nout,
            err_msg=f"data mismatch, {{ fout == nout }}. Failed to test {inspect.getsource(ex)}"
        )
