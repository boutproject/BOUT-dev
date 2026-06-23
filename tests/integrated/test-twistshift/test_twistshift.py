#!/usr/bin/env python3

import numpy as np
from boutdata import collect
from boututils.run_wrapper import launch_safe

datapath = "data"
nproc = 1
tol = 1.0e-13
non_periodic_tol = 1.0e-6

GUARD_PAIRS = [(0, -4), (1, -3), (2, -2), (3, -1)]


def test_twistshift():

    cmd = "./test-twistshift"

    s, out = launch_safe(cmd, nproc=nproc, pipe=True)

    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    test = collect("test", path=datapath, yguards=True, info=False)
    test_aligned = collect("test_aligned", path=datapath, yguards=True, info=False)
    result = collect("result", path=datapath, yguards=True, info=False)

    # from boututils.showdata import showdata
    # showdata([test, test_aligned, result], titles=['test', 'test_aligned', 'result'])

    # Check test_aligned is *not* periodic in y
    failures = []
    for ylower, yupper in GUARD_PAIRS:
        diff = np.abs(test_aligned[:, yupper, :] - test_aligned[:, ylower, :])
        if np.any(diff < non_periodic_tol):
            failures.append(
                f"Alignment Check: 'test_aligned' should not be periodic. "
                f"Slices at jy={yupper} and jy={ylower} are identical within {non_periodic_tol}."
            )

    # Check test and result are the same
    diff_result_test = np.abs(result - test)
    if np.any(diff_result_test > tol):
        max_diff = np.max(diff_result_test)
        failures.append(
            f"Fail - result has not been communicated correctly - is different from input. "
            f"Maximum divergence is {max_diff:.4e} (exceeds {tol:.4e} tolerance)."
        )

    # Check result is periodic in y
    for ylower, yupper in GUARD_PAIRS:
        diff = np.abs(result[:, yupper, :] - result[:, ylower, :])
        if np.any(diff > tol):
            max_diff = np.max(diff)
            failures.append(
                f"Fail - result should be periodic  jy={yupper} and jy={ylower} should not be different. "
                f"Maximum difference is {max_diff:.4e}."
            )

    assert not failures, "TwistShift validation encountered failures:\n" + "\n".join(
        failures
    )
