#!/usr/bin/env python3

#
# Run the test, check it completed successfully
#

from boututils.run_wrapper import launch_safe, shell

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


def test_invpar():
    flags = ""
    for i, f in enumerate(flags_src):
        fl = {"acoef": 1}
        fl.update(f)
        for k, v in fl.items():
            flags += f" {k}_{i}={v}"

    regions = ["", " mesh:ixseps1=0 mesh:ixseps2=0"]
    flags = [flags + r for r in regions]

    for nproc in [1, 2, 4]:
        cmd = "./test_invpar"

        print(f"   {nproc} processors....")
        for f in flags:
            shell(["rm data/BOUT.dmp.* 2> err.log"])

            # Run the case
            launch_safe(cmd + " -q -q -q " + f, nproc=nproc, mthread=1)
