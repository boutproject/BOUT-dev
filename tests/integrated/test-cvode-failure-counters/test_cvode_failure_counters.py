import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from boutdata.collect import collect
from boututils.datafile import DataFile

# CVODE diagnostics that are written with and without the option
CVODE_DIAGNOSTICS = [
    "cvode_nsteps",
    "cvode_nfevals",
    "cvode_nniters",
    "cvode_npevals",
    "cvode_nliters",
    "cvode_last_step",
    "cvode_last_order",
    "cvode_num_fails",
    "cvode_nonlin_fails",
    "cvode_stab_lims",
]

# Counters written with solver:cvode_failure_counters (cvode_precon_wtime, a
# wall time, is not compared)
COUNTERS = [
    "cvode_krylov_solves",
    "cvode_krylov_accepted",
    "cvode_krylov_fail_newton",
    "cvode_rhs_fail_newton",
    "cvode_newton_notconv",
    "cvode_step_solve_fails",
    "cvode_constraint_fails",
    "cvode_lin_conv_fails",
    "cvode_jtimes",
]


def run_case(assert_success_in_shell, name, *options):
    """Run the model with the input file in data/ and the given options, with
    output in the directory `name`"""
    path = Path(name)
    path.mkdir()
    shutil.copy(Path("data") / "BOUT.inp", path)
    assert_success_in_shell(
        f"./test_cvode_failure_counters -d {name} {' '.join(options)}"
    )
    return path


def read(path, names):
    return {name: collect(name, path=str(path), info=False) for name in names}


def variables(path):
    with DataFile(str(path / "BOUT.dmp.0.nc")) as datafile:
        return set(datafile.list())


@pytest.mark.parametrize("precon", ["none", "bbd"])
def test_option_does_not_change_solution(assert_success_in_shell, precon):
    off = run_case(
        assert_success_in_shell, "off", f"solver:cvode_precon_method={precon}"
    )
    on = run_case(
        assert_success_in_shell,
        "on",
        f"solver:cvode_precon_method={precon}",
        "solver:cvode_failure_counters=true",
    )

    # The solution, the number of RHS calls and the existing CVODE diagnostics
    # are bit-identical
    names = ["f", "t_array", "ncalls", "rhs_failures"] + CVODE_DIAGNOSTICS
    expected = read(off, names)
    actual = read(on, names)
    for name in names:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)

    # The counters are only written with the option
    assert not set(COUNTERS) & variables(off)
    assert set(COUNTERS) <= variables(on)


@pytest.mark.parametrize("precon", ["none", "bbd"])
def test_counters(assert_success_in_shell, precon):
    path = run_case(
        assert_success_in_shell,
        "counters",
        f"solver:cvode_precon_method={precon}",
        "solver:cvode_failure_counters=true",
    )
    c = read(path, COUNTERS + ["cvode_nonlin_fails", "rhs_failures"])

    # Every injected RHS failure stopped a Newton iteration
    assert c["rhs_failures"][-1] > 0
    np.testing.assert_array_equal(c["cvode_rhs_fail_newton"], c["rhs_failures"])

    # Every Krylov solve that did not converge was either accepted (reduced
    # residual in the first Newton iteration) or stopped a Newton iteration
    np.testing.assert_array_equal(
        c["cvode_lin_conv_fails"],
        c["cvode_krylov_accepted"] + c["cvode_krylov_fail_newton"],
    )
    assert np.all(c["cvode_krylov_solves"] >= c["cvode_lin_conv_fails"])

    # The other Newton failures are what remains of cvode_nonlin_fails
    assert np.all(c["cvode_newton_notconv"] >= 0)

    # No constraints are set, so every rejected step had a failed Newton solve,
    # and every failed Newton solve at least one failed Newton iteration
    np.testing.assert_array_equal(c["cvode_constraint_fails"], 0)
    assert c["cvode_step_solve_fails"][-1] > 0
    assert np.all(c["cvode_step_solve_fails"] <= c["cvode_nonlin_fails"])

    if precon == "none":
        # Without preconditioner the Krylov subspace (maxl = 2) is too small
        assert c["cvode_krylov_accepted"][-1] > 0
        assert c["cvode_krylov_fail_newton"][-1] > 0


def test_requires_newton_iteration():
    result = subprocess.run(
        "./test_cvode_failure_counters solver:cvode_failure_counters=true "
        "solver:func_iter=true",
        shell=True,
        check=False,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode != 0
    assert "requires Newton iteration" in result.stdout + result.stderr
