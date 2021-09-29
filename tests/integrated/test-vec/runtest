#!/usr/bin/env python3

# Cores: 4

from boututils.run_wrapper import build_and_log, launch_safe

build_and_log("vector communication test")

status, out = launch_safe("./testVec", nproc=4, pipe=True)

with open("run.log.4", "w") as f:
    f.write(out)

if status:
    print(" => Test failed")
    exit(status)

print(" => Test passed")
