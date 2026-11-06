#!/usr/bin/env python3
# Test `InvertPar`. Each BOUT++ run tests each case in `flags_src` below

import pytest
from boututils.run_wrapper import launch_safe

flags_src = [
    {"acoef": 1, "bcoef": 0, "ccoef": 0, "dcoef": 0, "ecoef": 0},
    {"acoef": 1.5, "bcoef": "'2.*sin(2*y)'"},
    {"acoef": 1},
    {"acoef": 1, "bcoef": 2},
    {"ccoef": 1.793},
    {"ccoef": 3, "bcoef": 0},
    {"dcoef": 3.5},
    {"ecoef": -1},
    {"input_field": "'ballooning(exp(-y*y)*cos(z)*gauss(x,0.2))'"},
]
flags = ""
# Construct variables with unique names
for i, f in enumerate(flags_src):
    fl = {"acoef": 1}
    fl.update(f)
    for k, v in fl.items():
        flags += f" {k}_{i}={v}"


@pytest.mark.parametrize("region", ["", "mesh:ixseps1=0 mesh:ixseps2=0"])
@pytest.mark.parametrize("nproc", [1, 2, 4])
def test_invpar(region, nproc):
    cmd = f"./test_invpar -q -q -q {flags} {region}"
    launch_safe(cmd, nproc=nproc, mthread=1)
