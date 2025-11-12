#!/usr/bin/env python3
import pytest
from boututils.run_wrapper import build_and_log, launch_safe
from boutdata.collect import collect

from numpy import max, abs

shift_types = ["shifted", "shiftedinterp"]
variables = [("ddy", "ddy_check"), ("ddy2", "ddy_check")]


@pytest.fixture(scope="module")
def build_project():
    """Build the project once per test module."""
    build_and_log("parallel slices test")


def run_case(shift_type, variable, build_project):
    s, out = launch_safe(
        "./test_yupdown mesh:paralleltransform:type=" + shift_type,
        nproc=1,
        pipe=True,
        verbose=True,
    )

    with open("run.log", "w") as f:
        f.write(out)

    v, v_check = variable
    print("Testing %s and %s ... " % (v, v_check))
    ddy = collect(v, path="data", xguards=False, yguards=False, info=False)
    ddy_check = collect(
        v_check, path="data", xguards=False, yguards=False, info=False
    )

    diff = max(abs(ddy - ddy_check))

    print(f"shifttype: {shift_type} Max difference {diff}")
    if diff < 2e-5:
        return True
    return False


@pytest.mark.parametrize("shift_type, variable", list(zip(shift_types, variables)))
def test_case(shift_type, variable, build_project):
    success = run_case(shift_type, variable, build_project)
    assert success, f"Test failed for shift_type={shift_type}: , variable={variable}"
