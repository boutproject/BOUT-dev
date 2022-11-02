#!/usr/bin/env python3

#
# Run the test, check it completed successfully
#

# Requires: netcdf
# Cores: 4

from __future__ import print_function

try:
    from builtins import str
except:
    pass
from boututils.run_wrapper import build_and_log, shell, launch
from boutdata.collect import collect
from sys import stdout, exit


build_and_log("Cyclic Reduction test")

flags = ["", "nsys=2", "nsys=5 periodic", "nsys=7 n=10"]

code = 0  # Return code
for nproc in [1, 2, 4]:
    cmd = "./test_cyclic"

    print("   %d processors...." % (nproc))
    r = 0
    for f in flags:
        stdout.write("\tflags '" + f + "' ... ")

        shell("rm data/BOUT.dmp.* 2> err.log")

        # Run the case
        status, out = launch(cmd + " " + f, nproc=nproc, mthread=1, pipe=True)
        with open(f"run.log.{nproc}.{r}", "w") as f:
            f.write(out)

        r = r + 1

        # Find out if it worked
        if status:
            print("PASSED")
        else:
            print("FAILED")
            code = 1

if code == 0:
    print(" => All cyclic reduction tests passed")
else:
    print(" => Some failed tests")

exit(code)
