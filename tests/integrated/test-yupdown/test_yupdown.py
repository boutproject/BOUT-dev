#!/usr/bin/env python3
import os
import pytest
from boututils.run_wrapper import launch_safe
from boutdata.collect import collect
from numpy import max, abs

shift_types = ["shifted", "shiftedinterp"]


@pytest.mark.parametrize("shift_type", shift_types)
def test_case(shift_type):

    os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "1"  # MPI oversubscribe

    s, out = launch_safe(
        "./test_yupdown mesh:paralleltransform:type=" + shift_type,
        nproc=1,
        pipe=True,
        verbose=True,
    )

    with open("run.log", "w") as f:
        f.write(out)

    errors = []
    for v, v_check in [("ddy", "ddy_check"), ("ddy2", "ddy_check")]:
        print("Testing %s and %s ... " % (v, v_check))
        ddy = collect(v, path="data", xguards=False, yguards=False, info=False)
        ddy_check = collect(
            v_check, path="data", xguards=False, yguards=False, info=False
        )

        diff = max(abs(ddy - ddy_check))

        print(f"shifttype: {shift_type} | {v} Max difference {diff}")
        if diff >= 2e-5:
            errors.append(f"{v} exceeded tolerance (diff: {diff})")

    assert not errors, f"Failures for {shift_type}:\n" + "\n".join(errors)
