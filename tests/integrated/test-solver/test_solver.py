#!/usr/bin/env python3

from boututils.run_wrapper import build_and_log, launch_safe

from sys import exit

nthreads = 1
nproc = 1


build_and_log("Solver test")

print("Running solver test")
status, out = launch_safe("./test_solver", nproc=nproc, mthread=nthreads, pipe=True)
with open("run.log", "w") as f:
    f.write(out)

if status:
    print(out)

exit(status)
