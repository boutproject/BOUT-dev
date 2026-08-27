import subprocess
from pathlib import Path

import numpy as np
import pytest
from boutdata.collect import collect

TEST_DIR = Path(".")
TIMESTEP = 0.2
CONSTRAINT_FACTOR = 0.5
EXPECTED_U = 1.0 / (1.0 + CONSTRAINT_FACTOR * TIMESTEP)
EXPECTED_PHI = CONSTRAINT_FACTOR * EXPECTED_U


def run_completed_with_known_mpi_finalize_issue(
    result: subprocess.CompletedProcess,
) -> bool:
    return (
        result.returncode == 143
        and "Run finished at" in result.stdout
        and "internal_Finalize" in result.stderr
        and "OFI poll failed" in result.stderr
    )


def run_command(command: str) -> str:
    result = subprocess.run(
        command, shell=True, capture_output=True, text=True, timeout=600, check=False
    )
    assert result.returncode == 0 or run_completed_with_known_mpi_finalize_issue(
        result
    ), f"Failed in {TEST_DIR}\nStderr: {result.stderr}\nOutput: {result.stdout}"
    return result.stdout + result.stderr


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
    run_command(f"./test_snes_constraints solver:equation_form={equation_form}")
    assert_constraint_solution()


def test_constraint_fieldsplit():
    output = run_command(
        "./test_snes_constraints "
        "solver:equation_form=backward_euler "
        "solver:pc_type=fieldsplit "
        "petsc:pc_fieldsplit_type=additive "
        "petsc:fieldsplit_diff_ksp_type=preonly "
        "petsc:fieldsplit_diff_pc_type=jacobi "
        "petsc:fieldsplit_alg_ksp_type=preonly "
        "petsc:fieldsplit_alg_pc_type=jacobi"
    )

    assert "Using PCFieldSplit preconditioner for DAE system" in output
    assert_constraint_solution()
