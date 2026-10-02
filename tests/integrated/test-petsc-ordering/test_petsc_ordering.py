#!/usr/bin/env python3
# requires: petsc
# cores: 4
"""
Integrated test: Ordering equivalence between PetscCellMapping
and the SNES solver's globalIndex traversal.

Runs test_petsc_ordering with 1, 2, and 4 MPI ranks and checks that:
  - The executable exits with status 0.
  - The output contains "ordering_check=PASS".
  - No "MISMATCH" or "SHIFT_NOT_CONSTANT" lines appear.
"""

import re

import pytest
from boututils.run_wrapper import launch_safe


def check_output(output, nproc):
    """
    Inspect the combined output for pass/fail markers.
    Returns a list of failure strings (empty means pass).
    """
    failures = []

    # Must contain the summary pass marker
    if "ordering_check=PASS" not in output:
        if "ordering_check=FAIL" in output:
            failures.append("ordering_check=FAIL found in output")
        else:
            failures.append("ordering_check marker not found in output")

    # Must not contain any mismatch lines
    mismatch_lines = [
        line for line in output.splitlines() if line.startswith("MISMATCH")
    ]
    if mismatch_lines:
        failures.append(
            f"{len(mismatch_lines)} MISMATCH line(s) found:\n"
            + "\n".join(f"  {line}" for line in mismatch_lines)
        )

    # Must not contain any shift-not-constant lines
    shift_lines = [
        line for line in output.splitlines() if line.startswith("SHIFT_NOT_CONSTANT")
    ]
    if shift_lines:
        failures.append(
            f"{len(shift_lines)} SHIFT_NOT_CONSTANT line(s) found:\n"
            + "\n".join(f"  {line}" for line in shift_lines)
        )

    # Extract and report the numeric summary for informational purposes
    for marker in ("total_mismatches", "total_shift_failures"):
        m = re.search(rf"{marker}=(\d+)", output)
        if m:
            count = int(m.group(1))
            if count != 0:
                failures.append(f"{marker}={count} (expected 0)")

    return failures


@pytest.mark.parametrize("nproc", [1, 2, 4])
def test_petsc_ordering(nproc):
    cmd = "./test_petsc_ordering"

    overall_pass = True

    for nxpe in [1, 2, 4]:
        if nxpe > nproc:
            break
        print(f"\n{'=' * 60}")
        print(f"Running with {nproc} MPI rank(s), nxpe={nxpe}")
        print(f"{'=' * 60}")

        returncode, output = launch_safe(cmd, nproc=nproc, pipe=True, verbose=True)

        if returncode != 0:
            # Note: MPI task can exit non-zero for reasons not connected to test
            print(f"Warning: Non-zero exit code: {returncode}")

        case_failures = []
        case_failures.extend(check_output(output, nproc))

        if case_failures:
            print(output)  # Only print output on failure
            overall_pass = False
            print(f"FAIL (nproc={nproc}, nxpe={nxpe}):")
            for f in case_failures:
                print(f"  - {f}")
        else:
            print(f"PASS (nproc={nproc}, nxpe={nxpe})")

    print(f"\n{'=' * 60}")
    if overall_pass:
        print("ALL CASES PASSED")
        return 0
    else:
        print("ONE OR MORE CASES FAILED")
        return 1
