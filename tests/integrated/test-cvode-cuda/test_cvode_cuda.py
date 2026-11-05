"""Requires a CUDA GPU; run on a GPU allocation, not a login node."""

import pytest
from boututils.run_wrapper import launch_safe


@pytest.mark.parametrize("nproc", [1, 2])
@pytest.mark.parametrize("backend", ["sundials", "cuda"])
def test_analytic_decay(backend, nproc):
    status, output = launch_safe(
        f"./test_cvode_cuda solver:nvector={backend}", nproc=nproc, pipe=True
    )
    assert status == 0, output
    if backend == "cuda":
        assert "SUNDIALS CUDA + MPIPlusX" in output
