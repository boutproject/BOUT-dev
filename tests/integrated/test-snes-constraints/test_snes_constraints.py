from pathlib import Path

import numpy as np
import pytest
from boutdata.collect import collect
from boututils.run_wrapper import launch

TEST_DIR = Path(".")
TIMESTEP = 0.2
CONSTRAINT_FACTOR = 0.5
# For one backward-Euler step:
#   (u1 - u0) / dt = -u1 + phi1,  phi1 = c * u1
# so u1 = u0 / (1 + (1 - c) * dt)
EXPECTED_U = 1.0 / (1.0 + (1.0 - CONSTRAINT_FACTOR) * TIMESTEP)
EXPECTED_PHI = CONSTRAINT_FACTOR * EXPECTED_U


def assert_constraint_solution():
    u = np.asarray(collect("u", path="data", info=False))[-1]
    phi = np.asarray(collect("phi", path="data", info=False))[-1]

    np.testing.assert_allclose(u, EXPECTED_U, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(phi, EXPECTED_PHI, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(phi - CONSTRAINT_FACTOR * u, 0.0, atol=1e-12, rtol=0.0)


@pytest.mark.parametrize(
    "equation_form",
    ["backward_euler", "rearranged_backward_euler", "pseudo_transient"],
)
def test_constraint_equation_forms(equation_form):
    status, output = launch(
        f"./test_snes_constraints solver:equation_form={equation_form}",
        nproc=1,
        pipe=True,
    )

    if status:
        print(f"WARNING: status = {status}")
        print(output)

    assert_constraint_solution()


def test_constraint_fieldsplit():
    status, output = launch(
        "./test_snes_constraints "
        "solver:equation_form=backward_euler "
        "solver:pc_type=fieldsplit "
        "petsc:pc_fieldsplit_type=additive "
        "petsc:fieldsplit_diff_ksp_type=preonly "
        "petsc:fieldsplit_diff_pc_type=jacobi "
        "petsc:fieldsplit_alg_ksp_type=preonly "
        "petsc:fieldsplit_alg_pc_type=jacobi",
        nproc=1,
        pipe=True,
    )

    if status:
        print(f"WARNING: status = {status}")
        print(output)

    assert "Using PCFieldSplit preconditioner for DAE system" in output
    assert_constraint_solution()
