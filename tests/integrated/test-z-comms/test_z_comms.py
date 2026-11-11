#!/usr/bin/env python3
#
# Run the Z communication test for several decompositions in X and Z

import numpy as np
from boututils.datafile import DataFile
from boututils.run_wrapper import launch_safe

# (number of processors, NXPE, NZPE, MZG)
CASES = [
    (1, 1, 1, 0),  # No Z guards: periodicity by index wrapping, no Z communication
    (1, 1, 1, 2),  # One processor in Z with guards: sends to itself
    (2, 1, 2, 2),
    (3, 1, 3, 2),
    (4, 1, 4, 2),
    (4, 2, 2, 2),  # Internal X boundaries: X-Z corners from X and Z communication
    (6, 2, 3, 3),
]


GRIDFILE = "z_comms_grid.nc"


def write_grid():
    """Grid file with the same sizes as data/BOUT.inp. The 3D variable has no Z
    guard cells, and values 1000 * x + z"""
    nx, ny, nz = 8, 1, 12
    x = np.arange(nx).reshape((nx, 1, 1))
    z = np.arange(nz).reshape((1, 1, nz))
    with DataFile(GRIDFILE, create=True, format="NETCDF4") as f:
        f.write("nx", nx)
        f.write("ny", ny)
        f.write("nz", nz)
        f.write("f3d", (1000.0 * x + z) * np.ones((nx, ny, nz)))


def test_z_comms():
    write_grid()
    failures = []
    for (nproc, nxpe, nzpe, mzg), gridfile in [
        (case, g) for case in CASES for g in (False, True)
    ]:
        cmd = f"./test_z_comms NXPE={nxpe} NZPE={nzpe} MZG={mzg}"
        if gridfile:
            cmd += f" mesh:file={GRIDFILE} check_gridfile=true"
        print(f"Running command: {cmd} on {nproc} processors")
        try:
            status, out = launch_safe(cmd, nproc=nproc, pipe=True)
        except RuntimeError as e:
            status, out = 1, str(e)
        with open(f"run.log.{nproc}.{nxpe}.{nzpe}.{mzg}.{int(gridfile)}", "w") as f:
            f.write(out)
        if status != 0 or "Z communication test passed" not in out:
            failures.append(f"{nproc=} {nxpe=} {nzpe=} {mzg=} {gridfile=}:\n{out}")

    assert not failures, "\n".join(failures)


if __name__ == "__main__":
    test_z_comms()
