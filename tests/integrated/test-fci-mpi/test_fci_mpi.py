#!/usr/bin/env python3
#
# Python script to run and analyse MPI test

import numpy.testing as npt
import pytest
from boutdata.collect import collect, create_cache
from boututils.run_wrapper import launch_safe

MAXCORES = 8
NSLICES = [1]
IMPLEMENTATIONS = ["hermitespline", "monotonichermitespline"]
COLLECT_KW = {"info": False, "xguards": False, "yguards": False, "path": "data"}


def run_case(nxpe: int, nype: int, mthread: int, nslice: int, implementation: str):
    cmd = f"./fci_mpi NXPE={nxpe} NYPE={nype} mesh:paralleltransform:xzinterpolation:type={implementation}"
    print(f"Running command: {cmd}")

    _, out = launch_safe(cmd, nproc=nxpe * nype, mthread=mthread, pipe=True)

    # Save output to log file
    with open(f"run.log.{nxpe}.{nype}.{nslice}.log", "w") as f:
        f.write(out)


@pytest.mark.parametrize("implementation", IMPLEMENTATIONS)
@pytest.mark.parametrize("nslice", NSLICES)
def test_fci_mpi(implementation, nslice, subtests):
    run_case(1, 1, MAXCORES, nslice, implementation)
    cache = create_cache(path="data", prefix="BOUT.dmp")

    ref = {}
    for i in range(4):
        for yp in range(1, nslice + 1):
            for y in [-yp, yp]:
                name = f"output_{i}_{y:+d}"
                ref[name] = collect(name, datafile_cache=cache, **COLLECT_KW)

    # Explicitly delete cache to close netCDF files
    del cache

    for nxpe, nype in [(2, 1), (1, 2), (2, 2)]:
        mthread = MAXCORES // (nxpe * nype)
        run_case(nxpe, nype, mthread, nslice, implementation)

        cache = create_cache(path="data", prefix="BOUT.dmp")
        results = {
            name: collect(name, datafile_cache=cache, **COLLECT_KW) for name in ref
        }
        # Explicitly delete cache to close netCDF files
        del cache

        for name, val in ref.items():
            with subtests.test(msg=f"{name}, {nxpe=}, {nype=}"):
                npt.assert_allclose(val, results[name])
