#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

import pathlib
import pytest
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect, create_cache

errors = [
    "max_error1",
    "max_error2",
    "max_error3",
    "max_error5",
    "max_error6",
    "max_error7",
]
tol = 2e-4  # Absolute (?) tolerance


@pytest.mark.parametrize("nproc, nzpe", [(1, 1), (2, 1), (2, 2), (4, 1)])
def test_petsc_laplace(nproc, nzpe):
    cmd = f"./test_petsc_laplace NZPE={nzpe}"

    shell(["rm data/BOUT.dmp.*.nc"])

    print(f"   {nproc} processors, {nzpe=}....")
    s, out = launch_safe(cmd, nproc=nproc, pipe=True, verbose=True)

    pathlib.Path(f"run.log.{nproc}").write_text(out)
    cache = create_cache(path="data", prefix="BOUT.dmp")

    # Collect output data
    for varname in errors:
        error = collect(varname, path="data", info=False, datafile_cache=cache)

        assert error > 0, f"Convergence error for {varname} on {nproc} proc(s)"
        assert error <= tol, (
            f"Fail: {varname} maximum error is {error:e} (Tolerance: {tol:e})"
        )
