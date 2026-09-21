#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

import pytest
from boututils.run_wrapper import shell, launch, launch_safe
from boutdata.collect import collect

tol_orth = 5.0e-8

# Note accuracy of test is limited when g12!=0 by inconsistency between the way boundary
# conditions are applied in LaplaceXY and the way they are applied in the D2DXDY()
# operator called by Laplace_perp(). In D2DXDY(f) 'free_o3' boundary conditions are
# applied to dfdy before calculating DDX(dfdy).
tol_nonorth = 2.0e-5

argslist = [
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=dirichlet "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=dirichlet f:bndry_ydown=dirichlet",
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=neumann "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann",
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=free_o3 "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=free_o3 f:bndry_ydown=free_o3",
    #'laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=dirichlet '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=dirichlet f:bndry_ydown=dirichlet',
    #'laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=neumann '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann',
    #'laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=dirichlet '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=dirichlet f:bndry_ydown=dirichlet',
    #'laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=neumann '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann',
    #'laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=dirichlet '
    #'f:bndry_xin=neumann f:bndry_xout=neumann f:bndry_yup=dirichlet f:bndry_ydown=dirichlet',
    "laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=neumann "
    "f:bndry_xin=neumann f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann laplacexy:pctype=hypre",
    "laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=free_o3 "
    "f:bndry_xin=neumann f:bndry_xout=dirichlet f:bndry_yup=free_o3 f:bndry_ydown=free_o3 laplacexy:pctype=hypre",
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=neumann "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann b:function=.1",
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=free_o3 "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=free_o3 f:bndry_ydown=free_o3 b:function=.1",
    "laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=neumann "
    "f:bndry_xin=neumann f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann b:function=.1 laplacexy:pctype=hypre",
    "laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=free_o3 "
    "f:bndry_xin=neumann f:bndry_xout=dirichlet f:bndry_yup=free_o3 f:bndry_ydown=free_o3 b:function=.1 laplacexy:pctype=hypre",
]


@pytest.mark.parametrize(
    "is_nonorth, tol",
    [(False, tol_orth), (True, tol_nonorth)],
    ids=["orthogonal", "non_orthogonal"],
)
@pytest.mark.parametrize("case_idx, base_args", enumerate(argslist))
def test_laplacexy(is_nonorth, tol, case_idx, base_args):
    nproc = 8

    args = base_args if is_nonorth else f"{base_args} mesh:g12=0."
    cmd = f"./test-laplacexy {args}"

    shell(["rm -f data/BOUT.dmp.*.nc"])

    if "hypre" in args:
        s, out = launch(cmd, nproc=nproc, pipe=True, verbose=True)
        if s == 134:
            # PETSc did not recognise pctype option, probably means it
            # was not compiled with hypre, so skip tests that need
            # hypre to converge
            pytest.skip("hypre not available as pre-conditioner in PETSc. Skipping...")
    else:
        s, out = launch_safe(cmd, nproc=nproc, pipe=True, verbose=True)

    label = "nonorth" if is_nonorth else "orth"
    log_filename = f"run.log.{nproc}.{label}.case_{case_idx}"
    with open(log_filename, "w") as f:
        f.write(out)

    # Collect output data
    error = collect("max_error", path="data", info=False)

    assert error > 0, "Convergence error"
    assert error <= tol, (
        f"Validation failure on case {case_idx} ({label}). "
        f"Maximum absolute error {error:.4e} exceeds defined tolerance threshold {tol:.4e}"
    )
