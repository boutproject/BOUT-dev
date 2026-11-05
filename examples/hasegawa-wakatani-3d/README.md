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

Profiling with NVIDIA NVTX
--------------------------

Rebuild BOUT++ and the example together to capture the full PVODE hierarchy:

```sh
cmake -S . -B build -DUSE_NVTX=ON -DBOUT_USE_PVODE=ON \
  -DNVTX3_INCLUDE_DIR=/path/to/nvtx/include
cmake --build build --target hw3d -j
cd build/examples/hasegawa-wakatani-3d
nsys profile --trace=nvtx,cuda,mpi -o hw3d ./hw3d -d data
```

Set `type = pvode` in the `[solver]` section of `data/BOUT.inp`.
Ranges follow `PvodeSolver::run` → `PvodeSolver::integrate` → `PVODE::CVode`
→ `CVStep` → `StepAttempt` → `CVnls` → `CVnlsNewton` → `CVNewtonIteration`
→ `NewtonIteration` → `CVSpgmrSolve` → `SpgmrSolve` → `GMRESIteration`.
Matrix-vector callbacks nest `PvodeSolver::rhs`, `Solver::run_rhs`, and the
HW3D model ranges beneath these levels. Colors cycle automatically with nesting
depth: blue, cyan, green, yellow, orange, magenta, purple, then blue again.
Depth is tracked separately for each host thread. Scopes close automatically on
returns and exceptions, and introduce no GPU synchronization.

`USE_NVTX` defaults to OFF. A standalone example built against an existing
uninstrumented BOUT++ installation records only the model ranges; rebuild the
library with this option for the solver levels.


CVODE CUDA vector operations
----------------------------

With BOUT++ CUDA and Umpire enabled, and a SUNDIALS installation exporting both
`SUNDIALS::nveccuda` and `SUNDIALS::nvecmpiplusx`, configure reports
`SUNDIALS CUDA/MPIPlusX vector backend enabled for CVODE`.
Reconfigure and rebuild on the GPU build host in the original module/Spack environment:

```sh
cmake -S . -B build_release_RAJA_CUDA_SUNGPU
cmake --build build_release_RAJA_CUDA_SUNGPU --target hw3d -j 8
```

In `data/BOUT.inp`, select the new backend explicitly:

```ini
[solver]
type = cvode
nvector = cuda
cvode_precon_method = none
```

The startup log must say
`Using SUNDIALS CUDA + MPIPlusX (GPU vector operations) N_Vector backend`.
`nvector=sundials` continues to select CPU parallel vectors. `nvector=cuda` is
currently implemented only for CVODE; PETSc/BBD preconditioners and Jacobian
export are rejected. User callbacks operate on BOUT++ fields as before.

CVODE's vector data live in CUDA device memory. GPU packing/unpacking kernels
connect them to BOUT++ fields with the same ordering and evolved boundary cells
as the CPU backend. BOUT++ field allocation is unchanged (currently Umpire pinned
memory); this does not move all model storage to GPU HBM. Packing synchronizes
with SUNDIALS and BOUT++ streams for correctness. Per-field tolerance and
constraint arrays are uploaded during setup; RHS state transfers do not stage
through a host solution vector. MPIPlusX performs inter-rank reductions without
requiring CUDA-aware MPI. Each rank must have its intended CUDA device selected
by the job environment (for example, `CUDA_VISIBLE_DEVICES`).

Profile a short run on a GPU allocation:

```sh
cd build_release_RAJA_CUDA_SUNGPU/examples/hasegawa-wakatani-3d
nsys profile --trace=cuda,nvtx -o hw3d-cvode-cuda \
  ./hw3d -d data solver:type=cvode solver:nvector=cuda solver:nout=1
nsys stats --report cuda_gpu_kern_sum hw3d-cvode-cuda.nsys-rep
```

Look for SUNDIALS CUDA vector arithmetic and reduction kernels in addition to
HW3D's RAJA and Laplacian kernels. `sundials_cuda_pack` kernels perform field/vector
transfers. NVTX host ranges alone do not demonstrate vector offload.

An analytic regression covers mixed 2D/3D fields, different evolved boundaries,
per-field tolerances, positivity constraints, and one/two MPI ranks. On a GPU
allocation, with tests enabled (`BOUT_TESTS=ON`), build target `test-cvode-cuda`, then run:

```sh
ctest --test-dir build_release_RAJA_CUDA_SUNGPU -R '^test-cvode-cuda$' --output-on-failure
```

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
