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
# Decompositions in Z: (NXPE, NYPE, NZPE, MZG). NZPE = 1 with MZG > 0 sends to itself
ZCASES = [(1, 1, 1, 2), (1, 1, 2, 2), (1, 1, 4, 2), (2, 2, 2, 2), (1, 4, 2, 2)]


def test_fci_mpi():

    COLLECT_KW = dict(info=False, xguards=False, yguards=False, path="data")

    def run_case(nxpe: int, nype: int, mthread: int, nzpe: int = 1, mzg: int = 0):

        cmd = f"./fci_mpi NXPE={nxpe} NYPE={nype} mesh:paralleltransform:xzinterpolation:type={implementation}"
        nproc = nxpe * nype
        if (nzpe, mzg) != (1, 0):
            cmd += f" NZPE={nzpe} MZG={mzg}"
            nproc *= nzpe
        print(f"Running command: {cmd}")

        _, out = launch_safe(cmd, nproc=nproc, mthread=mthread, pipe=True)

        # Save output to log file
        with open(f"run.log.{nxpe}.{nype}.{nzpe}.{mzg}.{nslice}.log", "w") as f:
            f.write(out)

    def test_case(
        nxpe: int, nype: int, mthread: int, ref: dict, nzpe: int = 1, mzg: int = 0
    ) -> bool:
        run_case(nxpe, nype, mthread, nzpe, mzg)

        failures = ""

        # Splitting in Z must not change the result beyond round-off
        rtol = 1e-7 if (nzpe, mzg) == (1, 0) else 1e-12

        for name, val in ref.items():
            try:
                npt.assert_allclose(val, collect(name, **COLLECT_KW), rtol=rtol)
            except AssertionError as e:
                failures += f"\ncase {nxpe=} {nype=} {nzpe=} {mzg=} {name=}\n{e}"

        return failures

    failures = ""

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
                failures += failures_

            if implementation == "monotonichermitespline":
                # Not implemented with guard cells in Z
                continue

            for nxpe, nype, nzpe, mzg in ZCASES:
                if nxpe * nype * nzpe > MAXCORES:
                    continue
                failures += test_case(nxpe, nype, 1, ref, nzpe, mzg)

        success = len(failures) == 0
        assert success, "\nSome tests failed:\n" + failures
