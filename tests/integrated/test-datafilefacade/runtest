#!/usr/bin/env python3
#
# Test the DataFileFacade interface by writing to dump files using the SAVE_ONCE macro
#
# requires: netcdf
# cores: 2

from boutdata.collect import collect
from boututils.run_wrapper import build_and_log, shell, launch_safe
import numpy
from sys import exit

build_and_log("DataFileFacade test")


success = True

testvars = {
    "f2d": 0.1,
    "f3d": 0.2,
    "fperp": 1.1,
    "v2d_contravariantx": 2.2,
    "v2d_contravarianty": 3.3,
    "v2d_contravariantz": 4.4,
    "v3d_contravariantx": 5.5,
    "v3d_contravarianty": 6.6,
    "v3d_contravariantz": 7.7,
    "v2d_covariant_x": 12.12,
    "v2d_covariant_y": 13.13,
    "v2d_covariant_z": 14.14,
    "v3d_covariant_x": 15.15,
    "v3d_covariant_y": 16.16,
    "v3d_covariant_z": 17.17,
    "integer": 42,
    "boolean": True,
    "real": 3.14,
    "name": "test string",
}


for nproc in [1, 2]:
    # delete any existing output
    shell("rm -f data/BOUT.dmp.*.nc data/BOUT.restart.*.nc")

    print(f"   {nproc} processor....")

    # run the test executable
    s, out = launch_safe("./test-datafile-facade", nproc=nproc, pipe=True)
    with open(f"run.log.{nproc}", "w") as f:
        f.write(out)

    # check the results
    for name, expected in testvars.items():
        # check non-evolving version
        result = collect(name, path="data", info=False)

        if result.dtype.kind in ("S", "U"):
            if str(result) != expected:
                success = False
                print(
                    f"{name} is different: got '{str(result)}', expected '{expected}'"
                )
            continue

        if not numpy.allclose(expected, result):
            success = False
            print(f"{name} is different: {numpy.max(numpy.abs(expected - result))}")

if success:
    print("=> All DataFileFacade tests passed")
    # clean up binary files
    shell(
        "rm -f data/BOUT.dmp.*.nc data/BOUT.restart.*.nc data/restart/BOUT.restart.0.nc"
    )
    exit(0)

print("=> Some failed tests")
exit(1)
