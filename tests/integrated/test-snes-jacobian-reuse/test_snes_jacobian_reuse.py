import re
from pathlib import Path

import numpy as np
import pytest
from boutdata.collect import collect
from boututils.run_wrapper import launch

TEST_DIR = Path(".")


def run_case(options):
    """Run the test model, returning the final state and the number of
    finite-difference Jacobian evaluations reported by PETSc."""
    log = TEST_DIR / "petsc_log.txt"
    log.unlink(missing_ok=True)

    status, output = launch(
        f"./test_snes_jacobian_reuse petsc:log_view=:{log} {options}",
        nproc=1,
        pipe=True,
    )
    if status:
        print(output)
    assert status == 0

    match = re.search(r"^MatFDColorApply\s+(\d+)", log.read_text(), re.MULTILINE)
    assert match is not None, "No finite-difference Jacobian evaluations in PETSc log"

    f = np.asarray(collect("f", path="data", info=False))[-1]
    g = np.asarray(collect("g", path="data", info=False))[-1]
    return f, g, int(match.group(1))


@pytest.mark.parametrize("matrix_free_operator", ["false", "true"])
def test_jacobian_reuse(matrix_free_operator):
    options = f"solver:matrix_free_operator={matrix_free_operator}"

    f_ref, g_ref, evals_ref = run_case(options)
    f, g, evals = run_case(options + " solver:jacobian_reuse=10")

    # Both runs approach the steady state f = g = 1
    for values in (f_ref, g_ref, f, g):
        np.testing.assert_allclose(values, 1.0, atol=1e-6, rtol=0.0)

    # Re-use must reduce the number of Jacobian evaluations substantially
    assert evals_ref >= 10
    assert 2 * evals <= evals_ref


@pytest.mark.parametrize(
    "options",
    [
        "solver:equation_form=pseudo_transient",
        "solver:equation_form=backward_euler",
        "solver:scale_rhs=true",
        "solver:use_coloring=false",
        "solver:matrix_free=true",
    ],
)
def test_jacobian_reuse_unsupported(options):
    """Combinations where the timestep term is not a simple diagonal shift
    are rejected rather than silently giving a wrong Jacobian."""
    status, output = launch(
        f"./test_snes_jacobian_reuse solver:jacobian_reuse=10 {options}",
        nproc=1,
        pipe=True,
    )
    assert status != 0
    assert "jacobian_reuse requires" in output
