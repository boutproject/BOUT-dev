"""Requires a CUDA GPU; run on a GPU allocation, not a login node."""

import pytest
from boututils.run_wrapper import launch_safe


@pytest.mark.parametrize("nproc", [1, 2])
@pytest.mark.parametrize("vector_fields", [False, True])
@pytest.mark.parametrize("monitor_timestep", [False, True])
@pytest.mark.parametrize(
    "backend,async_transfers",
    [("sundials", False), ("cuda", False), ("cuda", True)],
)
def test_analytic_decay(
    backend, async_transfers, nproc, vector_fields, monitor_timestep
):
    cuda_option = (
        f"solver:cuda_async_transfers={str(async_transfers).lower()} "
        if backend == "cuda"
        else ""
    )
    status, output = launch_safe(
        f"./test_cvode_cuda solver:nvector={backend} "
        f"{cuda_option}solver:monitor_timestep={str(monitor_timestep).lower()} "
        f"test:vector_fields={str(vector_fields).lower()}",
        nproc=nproc,
        pipe=True,
    )
    assert status == 0, output
    if backend == "cuda":
        assert "SUNDIALS CUDA + MPIPlusX" in output
        mode = "asynchronous" if async_transfers else "synchronous"
        assert f"CUDA callback transfers: {mode}" in output


@pytest.mark.parametrize("nproc", [1, 2])
@pytest.mark.parametrize("vector_fields", [False, True])
@pytest.mark.parametrize(
    "backend,async_transfers",
    [("sundials", False), ("cuda", False), ("cuda", True)],
)
def test_user_callbacks(backend, async_transfers, nproc, vector_fields):
    cuda_option = (
        f"solver:cuda_async_transfers={str(async_transfers).lower()} "
        if backend == "cuda"
        else ""
    )
    status, output = launch_safe(
        f"./test_cvode_cuda solver:nvector={backend} "
        f"{cuda_option}solver:cvode_precon_method=user solver:use_jacobian=true "
        f"test:user_callbacks=true test:vector_fields={str(vector_fields).lower()}",
        nproc=nproc,
        pipe=True,
    )
    assert status == 0, output
