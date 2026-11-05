#!/usr/bin/env python3
"""Compare resident GPU solves with CPU construction/solve on one GPU rank.

Usage: python check_resident_gpu.py /path/to/build/tests/integrated/test-laplace/test_laplace
Requires a CUDA or HIP BOUT++ build, numpy, and netCDF4. Includes varying A/C/D,
inhomogeneous boundary values, and gradient-boundary fallback via test_laplace.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import netCDF4
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--backend", choices=("cuda", "hip"), default="cuda")
    args = parser.parse_args()
    executable = args.executable.resolve()
    source = Path(__file__).resolve().parent
    option = "use_cusparse" if args.backend == "cuda" else "use_hipsparse"
    with tempfile.TemporaryDirectory(prefix="laplace-resident-") as tmp:
        outputs = []
        for gpu in (False, True):
            directory = Path(tmp) / ("gpu" if gpu else "cpu")
            directory.mkdir()
            (directory / "BOUT.inp").write_text((source / "data/BOUT.inp").read_text())
            command = [
                str(executable),
                "-d",
                str(directory),
                "NXPE=1",
                "laplace:type=cyclic",
                f"laplace:{option}={str(gpu).lower()}",
                f"mesh:file={source / 'test_laplace.grd.nc'}",
            ]
            with (directory / "run.log").open("w") as log:
                result = subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, check=False
                )
            if result.returncode:
                raise RuntimeError((directory / "run.log").read_text())
            with netCDF4.Dataset(directory / "BOUT.dmp.0.nc") as dataset:
                outputs.append(
                    {
                        name: np.asarray(value[:])
                        for name, value in dataset.variables.items()
                        if name.startswith("flag")
                    }
                )
        assert outputs[0].keys() == outputs[1].keys()
        for name in outputs[0]:
            np.testing.assert_allclose(
                outputs[1][name], outputs[0][name], rtol=1e-9, atol=1e-10, err_msg=name
            )
            print(
                f"{name}: max absolute difference "
                f"{np.max(np.abs(outputs[1][name] - outputs[0][name])):.3e}"
            )


if __name__ == "__main__":
    main()
