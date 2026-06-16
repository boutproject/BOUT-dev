#!/usr/bin/env python3
#
# Python script to run and analyse MPI test

from boututils.run_wrapper import launch_safe
from boutdata.collect import collect
import itertools

import numpy.testing as npt

# Resolution in x and y
NLIST = [1, 2, 4]
MAXCORES = 8
NSLICES = [1]


def test_fci_mpi():

    COLLECT_KW = dict(info=False, xguards=False, yguards=False, path="data")

    def run_case(nxpe: int, nype: int, mthread: int):

        cmd = f"./fci_mpi NXPE={nxpe} NYPE={nype} mesh:paralleltransform:xzinterpolation:type={implementation}"
        print(f"Running command: {cmd}")

        _, out = launch_safe(cmd, nproc=nxpe * nype, mthread=mthread, pipe=True)

        # Save output to log file
        with open(f"run.log.{nxpe}.{nype}.{nslice}.log", "w") as f:
            f.write(out)

    def test_case(nxpe: int, nype: int, mthread: int, ref: dict) -> bool:
        run_case(nxpe, nype, mthread)

        failures = []

        for name, val in ref.items():
            try:
                npt.assert_allclose(val, collect(name, **COLLECT_KW))
            except AssertionError as e:
                failures.append((nxpe, nype, name, e))

        return failures

    failures = []

    for implementation in ["hermitespline", "monotonichermitespline"]:
        for nslice in NSLICES:
            # reference data!
            run_case(1, 1, MAXCORES)

            ref = {}
            for i in range(4):
                for yp in range(1, nslice + 1):
                    for y in [-yp, yp]:
                        name = f"output_{i}_{y:+d}"
                        ref[name] = collect(name, **COLLECT_KW)

            for nxpe, nype in itertools.product(NLIST, NLIST):
                if (nxpe, nype) == (1, 1):
                    # reference case, done above
                    continue

                if nxpe * nype > MAXCORES:
                    continue

                mthread = MAXCORES // (nxpe * nype)
                failures_ = test_case(nxpe, nype, mthread, ref)
                failures.extend(failures_)

        success = len(failures) == 0

        assert success, "\nSome tests failed:"

        if not success:
            for nxpe, nype, name, error in failures:
                print("----------")
                print(f"case {nxpe=} {nype=} {name=}\n{error}")
