#!/usr/bin/env python3
#
# Run the Z communication test for several decompositions in X and Z

from boututils.run_wrapper import launch_safe

# (number of processors, NXPE, NZPE, MZG)
CASES = [
    (1, 1, 1, 2),  # One processor in Z with guards: sends to itself
    (2, 1, 2, 2),
    (3, 1, 3, 2),
    (4, 1, 4, 2),
    (4, 2, 2, 2),  # Internal X boundaries: X-Z corners from X and Z communication
    (6, 2, 3, 3),
]


def test_z_comms():
    failures = []
    for nproc, nxpe, nzpe, mzg in CASES:
        cmd = f"./test_z_comms NXPE={nxpe} NZPE={nzpe} MZG={mzg}"
        print(f"Running command: {cmd} on {nproc} processors")
        try:
            status, out = launch_safe(cmd, nproc=nproc, pipe=True)
        except RuntimeError as e:
            status, out = 1, str(e)
        with open(f"run.log.{nproc}.{nxpe}.{nzpe}.{mzg}", "w") as f:
            f.write(out)
        if status != 0 or "Z communication test passed" not in out:
            failures.append(f"{nproc=} {nxpe=} {nzpe=} {mzg=}:\n{out}")

    assert not failures, "\n".join(failures)


if __name__ == "__main__":
    test_z_comms()
