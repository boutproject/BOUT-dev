# BOUT++ CUDA Docker Container

This directory contains a Dockerfile for building BOUT++ with CUDA support.

## Changes from Previous Version

This Dockerfile addresses issue [#3233](https://github.com/boutproject/BOUT-dev/issues/3233) by updating dependencies:

- **Ubuntu**: 20.04 → 22.04 (better toolchain and CUDA support)
- **CUDA**: 12.2.2 → 12.6.1 (latest stable release)
- **Spack**: v0.21.1 → v0.23.0 (improved package management)
- **fmt**: 10.0.0 → 11.0.2 (fixes missing `fmt/base.h` header)
- **RAJA**: 2022.10.4 → 2024.07.0 (improved CUDA support)
- **Umpire**: 2022.03.1 → 2024.07.0 (improved memory management)
- **BLT**: 0.5.3 → 0.6.2 (build system updates)

## Prerequisites

- Docker or Podman installed
- NVIDIA Docker runtime (nvidia-docker2) for GPU access
- NVIDIA GPU with Compute Capability 8.0+ (Ampere or newer)
  - For other architectures, modify `cuda_arch` in the Dockerfile

## Building the Container

### Basic build:
```bash
cd .docker/cuda
docker build -t bout-cuda:latest .
```

### Build with specific BOUT++ commit:
```bash
docker build \
  --build-arg BOUT_URL=https://github.com/boutproject/BOUT-dev.git \
  --build-arg BOUT_COMMIT=<commit-hash-or-branch> \
  -t bout-cuda:custom .
```

## Running the Container

### Interactive shell:
```bash
docker run --gpus all -it bout-cuda:latest
```

### Run with mounted directory:
```bash
docker run --gpus all -it \
  -v $(pwd):/workspace \
  -w /workspace \
  bout-cuda:latest
```

### Check CUDA availability:
```bash
docker run --gpus all -it bout-cuda:latest bash -c "nvidia-smi && nvcc --version"
```

## GPU Architecture Notes

The Dockerfile is configured for NVIDIA A100 GPUs (Ampere, `cuda_arch=80`).

For other GPU architectures, modify these lines in the Dockerfile:

```dockerfile
# For Volta (V100): cuda_arch=70
spack add 'umpire@2024.07.0+cuda~examples+numa+openmp cuda_arch=70'
spack add 'raja@2024.07.0+cuda~examples~exercises+openmp cuda_arch=70'
# And in CMake configuration:
-DCUDA_ARCH="compute_70,code=sm_70"

# For Hopper (H100): cuda_arch=90
spack add 'umpire@2024.07.0+cuda~examples+numa+openmp cuda_arch=90'
spack add 'raja@2024.07.0+cuda~examples~exercises+openmp cuda_arch=90'
# And in CMake configuration:
-DCUDA_ARCH="compute_90,code=sm_90"
```

Common CUDA architectures:
- **70**: Volta (V100)
- **75**: Turing (RTX 20 series, T4)
- **80**: Ampere (A100, RTX 30 series)
- **86**: Ampere (RTX 30 mobile, A10, A40)
- **89**: Ada Lovelace (RTX 40 series)
- **90**: Hopper (H100)

## Troubleshooting

### Spack build failures
If Spack fails to build packages, try increasing Docker memory allocation (minimum 8GB recommended).

### CUDA runtime errors
Ensure nvidia-docker2 is properly installed:
```bash
# Test NVIDIA runtime
docker run --rm --gpus all nvidia/cuda:12.6.1-base-ubuntu22.04 nvidia-smi
```

### fmt header issues
If you encounter `fmt/base.h` not found errors, ensure you're using fmt 11.0+ (this is configured by default).

## Container Layout

- BOUT++ source: `/home/boutuser/BOUT-dev`
- BOUT++ install: `/opt/bout++`
- Spack environment: `/spack-env`
- Helper script: `/usr/local/bin/bout-env.sh` (loads Spack environment)

## Environment Activation

The container's entrypoint automatically activates the Spack environment. If you need to manually activate it:

```bash
. /spack/share/spack/setup-env.sh
spack env activate /spack-env
spack load cmake fmt raja umpire netcdf-cxx4 fftw
```

Or simply use the helper script:
```bash
/usr/local/bin/bout-env.sh <your-command>
```
