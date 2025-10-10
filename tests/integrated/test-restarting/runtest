#!/usr/bin/env python3

from boututils.run_wrapper import build_and_log, shell, launch_safe
from boutdata.collect import collect
import numpy as np
from sys import stdout, exit


build_and_log("restart test")

# Run once for 10 timesteps
s, out = launch_safe("./test_restarting solver:nout=10", nproc=1, pipe=True)

# Read reference data
f3d_0 = collect("f3d", path="data", info=False)
f2d_0 = collect("f2d", path="data", info=False)

###########################################
# Run twice, restarting and appending

print("-> Testing restart append")

shell("rm -f data/BOUT.dmp.0.nc")
s, out = launch_safe("./test_restarting solver:nout=5", nproc=1, pipe=True)
s, out = launch_safe(
    "./test_restarting solver:nout=5 restart append", nproc=1, pipe=True
)

f3d_1 = collect("f3d", path="data", info=False)
f2d_1 = collect("f2d", path="data", info=False)

success = True
tolerance = 1e-10

if f3d_1.shape != f3d_0.shape:
    print("Fail: Field3D field has wrong shape")
    success = False
if f2d_1.shape != f2d_0.shape:
    print("Fail: Field2D field has wrong shape")
    success = False

if not np.allclose(f3d_1, f3d_0, atol=tolerance):
    print("Fail: Field3D values differ")
    success = False

if not np.allclose(f2d_1, f2d_0, atol=tolerance):
    print("Fail: Field2D values differ")
    success = False

###########################################
# Test restart

print("-> Testing restart")

shell("rm -f data/BOUT.dmp.0.nc")
s, out = launch_safe("./test_restarting solver:nout=5", nproc=1, pipe=True)
s, out = launch_safe("./test_restarting solver:nout=5 restart", nproc=1, pipe=True)

f3d_1 = collect("f3d", path="data", info=False)
f2d_1 = collect("f2d", path="data", info=False)

if f3d_1.shape[0] != 6:
    print("Fail: Field3D has wrong shape")
    success = False
if f2d_1.shape[0] != 6:
    print("Fail: Field2D has wrong shape")
    success = False

if not np.allclose(f3d_1, f3d_0[5:, :, :, :], atol=tolerance):
    print("Fail: Field3D values differ")
    success = False
if not np.allclose(f2d_1, f2d_0[5:, :, :], atol=tolerance):
    print("Fail: Field2D values differ")
    success = False

if not success:
    print("=> Some tests failed")
    exit(1)

print("=> Success")
exit(0)
