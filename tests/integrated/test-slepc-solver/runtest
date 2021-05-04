#!/usr/bin/env python3

# requires: slepc

from boutdata.collect import collect
from boututils.run_wrapper import build_and_log, launch_safe
from numpy import isclose

build_and_log("SLEPc eigen solver test")

print("Running SLEPc eigen solver test")
status, out = launch_safe("./test-slepc-solver", nproc=1, pipe=True, verbose=True)

with open("run.log", "w") as f:
    f.write(out)

eigenvalues = collect("t_array", path="data", info=False)

expected_eigenvalues = [0.0, 1.0]

if isclose(expected_eigenvalues, eigenvalues).all():
    print(" => SLEPc test passed")
    exit(0)
else:
    print(" => SLEPc test failed")
    print("    Eigenvalues:", eigenvalues)
    exit(1)
