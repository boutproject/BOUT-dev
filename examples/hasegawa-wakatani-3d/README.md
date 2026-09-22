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
cmake -S . -B build_release_RAJA_HIP -DHW3D_USE_ROCTX=ON \
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
