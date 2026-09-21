#!/usr/bin/env python3

#
# Run the test, compare results
#

import pytest
from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect
from numpy import abs, seterr


@pytest.mark.parametrize("nproc", [1, 2, 4])
def test_fieldgroupcomm(nproc):
    # Good chance we'll do 0.0/0.0, which generates a warning
    # Ignore this warning
    seterr(divide="ignore", invalid="ignore")

    var_correct = "fld1"
    vars_comp = ["fld2", "fld3"]
    exe_name = "test_fieldgroupcomm"
    tol = 1e-10  # Relative tolerance

    cmd = f"./{exe_name}"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    print(f"Running FieldGroup comm test on {nproc} processor(s)...")
    _, out = launch_safe(cmd, nproc=nproc, pipe=True)
    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # Analyse result
    # /"Correct" answer
    f1 = collect(var_correct, path="data", info=False)
    f1max = abs(f1).max()
    # /Two different fields which should be identical to correct
    err = []
    for v in vars_comp:
        tmp = collect(v, path="data", info=False)
        relative_error = abs(f1 - tmp).max() / f1max

        if relative_error > tol:
            err.append(f"Field '{v}' mismatch: relative error is {relative_error:.4e})")

    assert not err, f"Failures encountered on {nproc} processor(s):\n" + "\n".join(err)
