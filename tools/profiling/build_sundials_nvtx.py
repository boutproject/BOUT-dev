#!/usr/bin/env python3
"""Build instrumented SUNDIALS 7.9.0 and relink an existing CUDA/MPI BOUT++ build.

Run on the target GPU build host with its compiler/MPI environment loaded.
The existing system SUNDIALS installation is left untouched.
"""

import argparse
import subprocess
from pathlib import Path

from instrument_sundials import patch

# Commit underlying the official v7.9.0 release tag.
REVISION = "312fc0f3684f27209ca9dc9249d194436eb41a7a"


def read_cache(path):
    values = {}
    for line in path.read_text().splitlines():
        if line and not line.startswith(("#", "//")) and "=" in line:
            key, value = line.split("=", 1)
            values[key.split(":", 1)[0]] = value
    return values


def run(*args):
    print("+", " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), check=True)


def build(directory, jobs):
    directory = directory.resolve()
    cache = read_cache(directory / "CMakeCache.txt")
    # Validate tools before downloading or changing the existing build.
    required = (
        "CMAKE_C_COMPILER",
        "CMAKE_CXX_COMPILER",
        "CMAKE_CUDA_COMPILER",
        "MPI_C_COMPILER",
        "MPI_CXX_COMPILER",
    )
    for key in required:
        if not cache.get(key) or not Path(cache[key]).is_file():
            raise ValueError(
                f"{key} unavailable: {cache.get(key, 'unset')}; "
                "run on the target CUDA/MPI build host"
            )
    if cache.get("BUILD_SHARED_LIBS") != "ON":
        raise ValueError("the NVTX bridge requires this helper's shared-library build")
    if cache.get("USE_ROCTX") == "ON":
        raise ValueError("disable USE_ROCTX before enabling USE_NVTX")
    nvtx = Path(cache.get("NVTX3_INCLUDE_DIR", ""))
    if not (nvtx / "nvtx3/nvToolsExt.h").is_file():
        raise ValueError("NVTX3_INCLUDE_DIR does not contain nvtx3/nvToolsExt.h")
    source = directory / "_deps/sundials-nvtx-src"
    output = directory / "_deps/sundials-nvtx-build"
    install = directory / "_deps/sundials-nvtx-install"
    if not source.exists():
        run(
            "git",
            "clone",
            "--depth",
            "1",
            "--branch",
            "v7.9.0",
            "https://github.com/LLNL/sundials.git",
            source,
        )
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != REVISION:
        raise ValueError(f"unexpected SUNDIALS checkout revision: {revision}")
    patch(source)
    options = [f"-D{key}={cache[key]}" for key in required]
    if cache.get("CMAKE_CUDA_ARCHITECTURES"):
        options.append(
            f"-DCMAKE_CUDA_ARCHITECTURES={cache['CMAKE_CUDA_ARCHITECTURES']}"
        )
    run(
        "cmake",
        "-S",
        source,
        "-B",
        output,
        *options,
        f"-DCMAKE_INSTALL_PREFIX={install}",
        "-DCMAKE_INSTALL_LIBDIR=lib",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_SHARED_LIBS=ON",
        "-DBUILD_STATIC_LIBS=OFF",
        "-DSUNDIALS_ENABLE_MPI=ON",
        "-DSUNDIALS_ENABLE_CUDA=ON",
        "-DSUNDIALS_INDEX_SIZE=64",
        "-DSUNDIALS_PRECISION=double",
        "-DEXAMPLES_ENABLE_C=OFF",
        "-DEXAMPLES_ENABLE_CXX=OFF",
        "-DEXAMPLES_INSTALL=OFF",
    )
    run("cmake", "--build", output, "--parallel", jobs)
    run("cmake", "--install", output)
    run(
        "cmake",
        "-S",
        cache["CMAKE_HOME_DIRECTORY"],
        "-B",
        directory,
        "-DUSE_NVTX=ON",
        "-DBOUT_DOWNLOAD_SUNDIALS=OFF",
        "-DBOUT_USE_SUNDIALS=ON",
        f"-DSUNDIALS_ROOT={install}",
        f"-DSUNDIALS_DIR={install / 'lib/cmake/sundials'}",
    )
    run("cmake", "--build", directory, "--target", "hw3d", "--parallel", jobs)
    print(
        "Rebuild complete. Rerun nsys to generate a new report; existing reports are unchanged."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_directory", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    try:
        build(args.build_directory, args.jobs)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
