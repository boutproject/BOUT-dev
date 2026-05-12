#!/usr/bin/env python3

#
# Run the test, compare results
#

# Requires: all_tests
# Requires: netcdf
# Requires: not metric_3d
# Cores: 4


from boututils.run_wrapper import shell, launch_safe
from boutdata.collect import collect
from numpy import abs, seterr


def test_fieldgroupcomm():

    # Good chance we'll do 0.0/0.0, which generates a warning
    # Ignore this warning
    seterr(divide="ignore", invalid="ignore")

    varCorrect = "fld1"
    varsComp = ["fld2", "fld3"]
    name = "FieldGroup comm"
    exeName = "test_fieldgroupcomm"
    tol = 1e-10  # Relative tolerance

    print("Running {nm} test".format(nm=name))
    success = True

    for nproc in [1, 2, 4]:
        nxpe = 1
        if nproc > 2:
            nxpe = 2

        cmd = "./{exe} ".format(exe=exeName)

        shell("rm data/BOUT.dmp.*.nc")

        print("   %d processors ...." % (nproc))
        s, out = launch_safe(cmd, nproc=nproc, pipe=True)
        with open("run.log." + str(nproc), "w") as f:
            f.write(out)

        # Analyse result
        # /"Correct" answer
        f1 = collect(varCorrect, path="data", info=False)
        f1max = abs(f1).max()
        # /Two different fields which should be identical to correct
        err = []
        for v in varsComp:
            tmp = collect(v, path="data", info=False)
            err.append(abs((f1 - tmp)).max() / f1max)

        for i, e in enumerate(err):
            if e > tol:
                print(
                    "Fail, in {i}th comparison relative error is {re}".format(i=i, re=e)
                )
                success = False

    assert success, " => Some failed tests"
