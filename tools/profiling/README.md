# Comparing NVTX and ROCTX region timings

`extract_region_timings.py` compares every recorded region across two or more
profiles (one input also works). It uses Python 3 and the standard library.

```bash
python3 tools/profiling/extract_region_timings.py nsys-report.txt rocprof-run/ --labels NVIDIA AMD --calls
python3 tools/profiling/extract_region_timings.py first_marker_api_stats.csv second_nvtx_sum.csv --csv > comparison.csv
```

Inputs can be mixed:

- Nsight Systems `nvtx_sum` text reports, including reports embedded in job logs.
- Nsight Systems `nvtx_sum` CSV reports with `Range`, `Instances`, and
  `Total Time (ns)` columns.
- rocprofv3 `*_marker_api_stats.csv` files with `Name`, `Calls`, and
  `TotalDurationNs` columns.
- Run directories, searched recursively. The script first looks for rocprofv3
  marker statistics, then `*nvtx_sum*.csv`, then `.txt`, `.log`, or `HW3D-*`
  text reports. Only the first available format is used to avoid counting both
  text and CSV exports of the same run. Pass individual files when selecting
  reports from a directory containing several formats.

Each directory is **one comparison input**: keep separate runs in separate
inputs. Files from MPI ranks within a run are combined as total duration divided
by total calls, rather than averaging rank averages. Results are in milliseconds;
speedup is the first input's average divided by the current input's average.
Missing regions and undefined ratios are `N/A`. `--calls` shows sample counts;
`--labels` overrides the default file/directory names.

Names are matched exactly, except that Nsight's empty-domain leading `:` is
removed. Named NVTX domains are preserved. Nested ranges overlap: these are
inclusive host range timings, not exclusive times or GPU kernel timings.

Raw `.nsys-rep`, SQLite, and rocprofv3 trace/JSON files are not read directly;
export the region summaries above first.

Run the regression checks with:

```bash
python3 -m unittest discover -s tools/profiling -p 'test_*.py'
```

## SUNDIALS CVODE NVTX ranges

BOUT++ built with `USE_NVTX=ON` instruments `CvodeSolver::run`,
`CvodeSolver::integrate`, RHS, preconditioner, and Jacobian-vector callbacks.
For internal solver levels, SUNDIALS itself must also be rebuilt.

On the GH200 build host, with the original CUDA/MPI environment loaded, run:

```bash
python3 tools/profiling/build_sundials_nvtx.py build_release_RAJA_CUDA_SUNGPU --jobs 4
```

The helper downloads the pinned official SUNDIALS 7.9.0 release, applies the
patch, builds shared CUDA/MPI libraries using the compilers and GPU architecture
from the existing BOUT++ cache, and installs under that build's
`_deps/sundials-nvtx-install`. It then reconfigures BOUT++ to use that installation
and rebuilds `hw3d`. The shared installation under `programs/` is untouched.
This helper targets the current double-precision, 64-bit-index CUDA/MPI build.
It requires the original target-host compilers; it cannot rebuild GH200 binaries
on an x86 login host without that toolchain.

For a separately managed SUNDIALS source/build, apply only the patch with:

```bash
python3 tools/profiling/instrument_sundials.py /path/to/sundials-7.9.0-source
```

Rebuild and install SUNDIALS with its original options, then configure BOUT++
with `USE_NVTX=ON` and `SUNDIALS_DIR` pointing at the patched installation.
The patch supports version 7.9.0 and can be applied repeatedly. Its C scopes
require GCC/Clang cleanup attributes and ELF weak symbols (Linux). No SUNDIALS
built-in profiling option or NVTX dependency is required in SUNDIALS itself.
The optional bridge exported by BOUT++ supplies NVTX calls and the shared
thread-local depth counter. Without this bridge the SUNDIALS ranges are inert.

The resulting hierarchy includes:

```text
CvodeSolver::run
  CvodeSolver::integrate
    SUNDIALS::CVode
      SUNDIALS::CVStep
        SUNDIALS::StepAttempt
          SUNDIALS::CVNls
            SUNDIALS::NewtonSolve
              SUNDIALS::NewtonIteration
                SUNDIALS::CVLsSolve
                  SUNDIALS::SPGMRSolve
                    SUNDIALS::GMRESIteration
                      SUNDIALS::CVLsATimes
                        SUNDIALS::CVLsDQJtimes
                          CvodeSolver::rhs
                            Solver::run_rhs
                              HW3D::rhs
```

Actual nesting depends on the selected linear/nonlinear solver and callbacks;
these internal iteration scopes target Newton and SPGMR (GMRES).
Colors repeat by depth across all these scopes: blue, cyan, green, yellow,
orange, magenta, purple. Returns, retry `continue`s, and loop `break`s close
ranges automatically. No synchronization or solver algorithm changes are added.
At startup CVODE reports whether both wrapper and internal ranges are enabled.

If the job environment's `LD_LIBRARY_PATH` includes the original SUNDIALS
installation, prepend the new installation's `lib` directory before running:

```bash
# From the BOUT++ repository root:
export LD_LIBRARY_PATH="$PWD/build_release_RAJA_CUDA_SUNGPU/_deps/sundials-nvtx-install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

Rerun Nsight Systems after rebuilding; existing `.txt` and `.nsys-rep` files
retain their old contents. For example, from the example build directory:

```bash
nsys profile --trace=nvtx,cuda,mpi --stats=true -o profile-hw3d-nvtx ./hw3d -d data
```

### Vector arithmetic ranges

The patch also adds `SUNDIALS::N_V*` ranges at the generic N_Vector dispatch
layer for standard arithmetic, norms/dot products, fused operations, vector-array
operations, and local/single-buffer reductions. Examples include
`SUNDIALS::N_VLinearSum`, `SUNDIALS::N_VScale`, `SUNDIALS::N_VDotProd`,
`SUNDIALS::N_VWrmsNorm`, and `SUNDIALS::N_VLinearCombination`. The same names
cover CPU parallel vectors and CUDA/MPIPlusX vectors. Fused fallback work stays
inside the fused operation range. Vector lifecycle/accessor and buffer-transfer
functions are excluded.

These are inclusive host durations: CPU arithmetic is measured directly, while
CUDA operations may only enqueue kernels; reductions can include waiting and MPI
communication. No extra device synchronization is inserted. Use the CUDA timeline
(or Nsight's NVTX GPU projection reports) for associated GPU execution time.
Nested `N_V*` ranges can occur with composite vectors; do not sum parent and child
ranges. Compare outermost vector ranges within CVode, excluding setup, and report
RHS/iteration counts alongside totals. Rebuild and reprofile both backends.

### CPU validation of the SUNDIALS patch

For an isolated check without CUDA or MPI, configure the patched source with
`SUNDIALS_ENABLE_MPI=OFF`, `SUNDIALS_ENABLE_CUDA=OFF`, `BUILD_SHARED_LIBS=ON`,
`BUILD_STATIC_LIBS=OFF`, and `CMAKE_BUILD_TYPE=Release`, then build it. Run the
integration check below, replacing the two SUNDIALS paths. It records calls with
a test NVTX header, checks depth colors and all iteration scopes during a real
stiff solve, verifies the solution, and checks balanced scopes after a deliberate
`CV_TOO_MUCH_WORK` failure (which prints an expected SUNDIALS error).

```bash
SUNDIALS_SOURCE=/path/to/patched-source
SUNDIALS_BUILD=/path/to/validation-build
# Run from the BOUT++ repository root; keep assertions enabled in this check.
g++ -std=c++17 -DUSE_NVTX -rdynamic \
  -Iinclude -Itools/profiling/tests \
  -I"$SUNDIALS_SOURCE/include" -I"$SUNDIALS_BUILD/include" \
  tools/profiling/tests/sundials_nvtx_integration.cxx \
  -L"$SUNDIALS_BUILD/src/cvode" -Wl,-rpath,"$SUNDIALS_BUILD/src/cvode" \
  -L"$SUNDIALS_BUILD/src/nvector/serial" -Wl,-rpath,"$SUNDIALS_BUILD/src/nvector/serial" \
  -L"$SUNDIALS_BUILD/src/sunlinsol/spgmr" -Wl,-rpath,"$SUNDIALS_BUILD/src/sunlinsol/spgmr" \
  -L"$SUNDIALS_BUILD/src/sundials" -Wl,-rpath,"$SUNDIALS_BUILD/src/sundials" \
  -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsolspgmr -lsundials_core \
  -o sundials-nvtx-check
./sundials-nvtx-check
```
