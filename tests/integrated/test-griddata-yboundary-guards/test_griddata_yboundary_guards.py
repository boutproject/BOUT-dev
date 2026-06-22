#!/usr/bin/env python3

# Cores: 6
# Requires: netcdf

import pytest
import numpy
from pathlib import Path
from netCDF4 import Dataset
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect


@pytest.mark.parametrize("topology", ["doublenull", "singlenull"])
@pytest.mark.parametrize("n_yguards", [0, 1, 2])
def test_griddata_yboundary_guards(topology, n_yguards):
    nx = 4
    ny = 24
    blocksize = ny // 6
    nproc = 6


    datadir = Path(f"data-{topology}-{n_yguards}")
    datadir.mkdir(parents=True, exist_ok=True)

    gridname = f"grid-{topology}-{n_yguards}.nc"
    file_path = datadir / gridname

    # first generate grid file to test
    with Dataset(file_path, "w") as gridfile:
        gridfile.createDimension("x", nx)

        if topology == "doublenull":
            gridfile.createDimension("y", ny + 4 * n_yguards)
        else:
            gridfile.createDimension("y", ny + 2 * n_yguards)

        gridfile.createVariable("nx", numpy.int32)
        gridfile["nx"][...] = nx

        gridfile.createVariable("ny", numpy.int32)
        gridfile["ny"][...] = ny

        gridfile.createVariable("y_boundary_guards", numpy.int32)
        gridfile["y_boundary_guards"][...] = n_yguards

        gridfile.createVariable("MXG", numpy.int32)
        gridfile["MXG"][...] = 1

        gridfile.createVariable("MYG", numpy.int32)
        gridfile["MYG"][...] = 2 if n_yguards == 0 else n_yguards

        gridfile.createVariable("ixseps1", numpy.int32)
        gridfile["ixseps1"][...] = nx // 2 - 1

        gridfile.createVariable("ixseps2", numpy.int32)
        gridfile["ixseps2"][...] = nx // 2 - 1

        gridfile.createVariable("jyseps1_1", numpy.int32)
        gridfile["jyseps1_1"][...] = blocksize - 1

        if topology == "doublenull":
            gridfile.createVariable("jyseps2_1", numpy.int32)
            gridfile["jyseps2_1"][...] = 2 * blocksize - 1

            gridfile.createVariable("ny_inner", numpy.int32)
            gridfile["ny_inner"][...] = 3 * blocksize

            gridfile.createVariable("jyseps1_2", numpy.int32)
            gridfile["jyseps1_2"][...] = 4 * blocksize - 1
        else:  # singlenull
            gridfile.createVariable("jyseps2_1", numpy.int32)
            gridfile["jyseps2_1"][...] = ny // 2

            gridfile.createVariable("ny_inner", numpy.int32)
            gridfile["ny_inner"][...] = ny // 2

            gridfile.createVariable("jyseps1_2", numpy.int32)
            gridfile["jyseps1_2"][...] = ny // 2

        gridfile.createVariable("jyseps2_2", numpy.int32)
        gridfile["jyseps2_2"][...] = 5 * blocksize - 1

        y_dim = ny + (4 * n_yguards if topology == "doublenull" else 2 * n_yguards)
        testdata = numpy.zeros([nx, y_dim])
        testdata[:, :] = numpy.arange(y_dim)[numpy.newaxis, :]
        gridfile.createVariable("test", float, ("x", "y"))
        gridfile["test"][...] = testdata

    shell([f"rm -f {datadir}/BOUT.dmp.*.nc run.log.*"])

    print(
        f"Running {topology} test with n_yguards={n_yguards} on {nproc} processors..."
    )

    s, out = launch_safe(f"./test_griddata -d {datadir}", nproc=nproc, pipe=True)

    with open(f"run.log.{topology}.{n_yguards}.{nproc}", "w") as f:
        f.write(out)

    testfield = collect("test", path=str(datadir), info=False, yguards=True)

    if topology == "doublenull":
        if n_yguards == 0:
                # output has 2 y-guard cells, but grid file did not
            myg = 2
            checkfield = list(numpy.zeros(myg))
            checkfield += list(numpy.arange(ny // 2))
            checkfield += list(numpy.arange(ny // 2) + checkfield[-1] + 1)
            checkfield += list(numpy.zeros(myg) + checkfield[-1])
        else:
            checkfield = []
            checkfield += list(numpy.arange(n_yguards))
            checkfield += list(numpy.arange(ny // 2) + checkfield[-1] + 1)
            checkfield += list(
                numpy.arange(ny // 2) + checkfield[-1] + 1 + 2 * n_yguards
            )
            checkfield += list(numpy.arange(n_yguards) + checkfield[-1] + 1)
    else:  # singlenull
        if n_yguards == 0:
                # output has 2 y-guard cells, but grid file did not
            myg = 2
            checkfield = list(numpy.zeros(myg))
            checkfield += list(numpy.arange(ny))
            checkfield += list(numpy.zeros(myg) + checkfield[-1])
        else:
            checkfield = []
            checkfield += list(numpy.arange(n_yguards))
            checkfield += list(numpy.arange(ny) + checkfield[-1] + 1)
            checkfield += list(numpy.arange(n_yguards) + checkfield[-1] + 1)

    checkfield = numpy.array(checkfield)

    # Test value of testfield
    max_diff = numpy.max(numpy.abs(testfield - checkfield))
    assert max_diff <= 1e-13, (
        f"testfield mismatch in {topology} configuration for n_yguards={n_yguards}. "
        f"Maximum difference encountered: {max_diff:.4e}"
    )
