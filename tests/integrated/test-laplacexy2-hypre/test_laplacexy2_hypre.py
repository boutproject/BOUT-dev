#!/usr/bin/env python3

#
# Run the test, compare results against the benchmark
#

# requires: hypre
# cores: 8

import pytest
from boutdata.collect import collect
from boututils.run_wrapper import launch_safe, shell

tol = 5.0e-8
nproc = 8

argslist = [
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry_dirichlet=true "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=dirichlet f:bndry_ydown=dirichlet",
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=neumann "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann",
    #'laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry_dirichlet=true '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=dirichlet f:bndry_ydown=dirichlet',
    #'laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry_dirichlet=true '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann',
    #'laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry_dirichlet=true '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=dirichlet f:bndry_ydown=dirichlet',
    #'laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry_dirichlet=true '
    #'f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann',
    #'laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry_dirichlet=true '
    #'f:bndry_xin=neumann f:bndry_xout=neumann f:bndry_yup=dirichlet f:bndry_ydown=dirichlet',
    "laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=neumann "
    "f:bndry_xin=neumann f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann",
    "laplacexy:core_bndry_dirichlet=true laplacexy:pf_bndry_dirichlet=true laplacexy:y_bndry=neumann "
    "f:bndry_xin=dirichlet f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann b:function=.1",
    "laplacexy:core_bndry_dirichlet=false laplacexy:pf_bndry_dirichlet=false laplacexy:y_bndry=neumann "
    "f:bndry_xin=neumann f:bndry_xout=dirichlet f:bndry_yup=neumann f:bndry_ydown=neumann b:function=.1",
]


@pytest.mark.parametrize("args", argslist)
def test_laplacexy2_hypre(args):

    print("   %d processors...." % nproc)

    cmd = f"./test-laplacexy {args}"

    shell(["rm data/BOUT.dmp.*.nc > /dev/null 2>&1"])

    s, out = launch_safe(cmd, nproc=nproc, pipe=True, verbose=True)

    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # Collect output data
    error = collect("max_error", path="data", info=False)

    assert error > 0, f"Convergence error: maximum error is {error}"
    assert error <= tol, (
        f"Fail: maximum error is {error}, which exceeds tolerance {tol}"
    )
