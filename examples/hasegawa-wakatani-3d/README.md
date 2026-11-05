3D Hasegawa-Wakatani model
==========================

Adds a simple Ohm's law to the 2D Hasegawa-Wakatani model,
so that second derivatives along the magnetic field appear in the model.
A 3D grid is then used to evolve drift-wave turbulence.


Profiling with ROCm ROCTx
------------------------

ROCTx is ROCm's equivalent of NVTX. Enable the optional host-side ranges in
an existing BOUT++ build (with examples enabled):

```sh
cmake -S . -B build_release_RAJA_HIP -DUSE_ROCTX=ON \
  -Drocprofiler-sdk-roctx_ROOT=/path/to/rocm
cmake --build build_release_RAJA_HIP --target hw3d -j
cd build_release_RAJA_HIP/examples/hasegawa-wakatani-3d
rocprofv3 --marker-trace --hip-trace --kernel-trace --output-format pftrace -- ./hw3d -d data
```

Use your ROCm installation prefix in place of `/path/to/rocm`. The same option
works when configuring this example as a standalone project against an installed
BOUT++. It defaults to OFF, so ordinary builds need no ROCTx dependency.

The timeline contains `HW3D::init` and `HW3D::rhs`, with nested RHS ranges for
`solve_phi`, `phi_minus_n`, `communicate`, and `evaluate_rhs`. The latter includes
accessor construction and a nested `rhs_kernel` range for the RAJA loop.
Ranges are balanced on exceptions. No GPU synchronization is added: ranges
measure host scopes and kernel submission; use the GPU kernel trace to inspect
actual device execution times.

See the [ROCTx documentation](https://rocm.docs.amd.com/projects/rocprofiler-sdk/en/latest/how-to/using-rocprofiler-sdk-roctx.html)
for the annotation API.

Profiling the solver
--------------------

Configure and rebuild BOUT++ itself with `-DUSE_NVTX=ON` (or `-DUSE_ROCTX=ON`)
to trace work between model RHS calls. Enabling profiling only in a standalone
example does not instrument an already installed library.

PVODE and CVODE ranges enclose integration and RHS callbacks. Nested ranges
identify solution/derivative vector copies and `Solver::pre_rhs`/`post_rhs`
boundary processing. Time inside integration but outside callbacks belongs to
the integrator; inspect CUDA and CPU tracks to distinguish GPU work and waits.

NVTX ranges use a repeating depth palette from `nvtxColor`: Blue, Cyan, Green,
Yellow, Orange, Magenta, Purple. Sibling scopes use the same color. Existing
`HW3D::` region names are unchanged. ROCTX records the same scopes without custom
colors. Profiling adds no GPU synchronization.

PVODE internals are also instrumented: `PVODE::CVStep`/`StepAttempt` enclose
step attempts; `CVPredict`, `CVEwtSet`, `CVRescale`, `CVDoErrorTest`, and
`CVCompleteStep` identify prediction, weighting, error checks, and history
updates. Nonlinear and GMRES iteration ranges enclose their solver operations.
`PVODE::N_V*` ranges measure CPU vector operations, with a nested
`PVODE::PVecAllReduce` range for MPI reductions. These scopes use the same
depth-based colors as the model and solver wrappers.

GPU Laplacian inversion
-----------------------

The resident GPU FFT path also handles nonperiodic X with second-order Dirichlet
boundaries (`inner_boundary_flags`/`outer_boundary_flags` containing only
`INVERT_BNDRY_ONE`, `INVERT_SET`, or `INVERT_RHS`, or zero). X must have one MPI
partition. Other boundary operators retain CPU matrix construction. Coefficients and metrics
remain in device buffers; content comparisons avoid redundant uploads while
detecting changes to A/C/D and in-place metric edits. Mapped CUDA fields use GPU
packing and unpacking to avoid host staging of the RHS and solution. Both spectral
transposes and matrix construction stay on the GPU. The optional first-call
comparison checks against CPU matrix and boundary-RHS construction.

A single-GPU regression compares the existing Laplacian test's coefficient and
boundary cases against CPU construction/solve, including gradient fallback:

```sh
python tests/integrated/test-laplace/check_resident_gpu.py \
  build_release_RAJA_CUDA_SUNGPU/tests/integrated/test-laplace/test_laplace
```
