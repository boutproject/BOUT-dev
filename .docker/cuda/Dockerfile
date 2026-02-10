# BOUT++ CUDA-enabled Docker container
# This Dockerfile builds BOUT++ with CUDA support using Spack for dependency management
# Fixes issue #3233 by using updated fmt library version

FROM docker.io/nvidia/cuda:12.6.1-devel-ubuntu22.04 AS base
MAINTAINER BOUT++ Development Team
USER root

# Install system dependencies
# Using Ubuntu 22.04 for more recent toolchain and better CUDA support
RUN \
 apt update &&\
 DEBIAN_FRONTEND=noninteractive TZ=America/Los_Angeles apt install -y \
   build-essential \
   git \
   software-properties-common \
   python3 \
   python3-pip \
   wget \
   mpich \
   libmpich-dev \
   unzip \
   libnuma-dev \
   libncurses-dev \
   curl \
   libcurl4-openssl-dev \
   gfortran \
   ca-certificates &&\
 apt upgrade -y &&\
 apt clean

FROM base AS setup-spack

# Install Spack and set up environment
RUN \
 git clone --depth 1 --single-branch --branch v0.23.0 https://github.com/spack/spack.git &&\
 . /spack/share/spack/setup-env.sh &&\
 spack env create -d /spack-env &&\
 spack env activate /spack-env &&\
 spack add cmake@3.24: &&\
 spack add blt@0.6.2 &&\
 spack add 'umpire@2024.07.0+cuda~examples+numa+openmp cuda_arch=80' &&\
 spack add 'raja@2024.07.0+cuda~examples~exercises+openmp cuda_arch=80' &&\
 spack add fmt@11.0.2 &&\
 spack add netcdf-cxx4@4.3.1 &&\
 spack add fftw@3.3.10 &&\
 spack external find --not-buildable &&\
 spack external find perl cuda bzip2 pkgconfig ncurses mpich curl --not-buildable &&\
 spack install -j$(nproc)

# Set up environment variables for the build
ENV SPACK_ROOT=/spack
ENV PATH="${SPACK_ROOT}/bin:${PATH}"

FROM setup-spack AS bout-build

# Create boutuser (for consistency with existing Docker setup)
RUN useradd -m -s /bin/bash boutuser && \
    echo "boutuser:boutforever" | chpasswd && \
    usermod -aG sudo boutuser

USER boutuser
WORKDIR /home/boutuser

# Clone and build BOUT++
ARG BOUT_URL=https://github.com/boutproject/BOUT-dev.git
ARG BOUT_COMMIT=main

RUN git clone ${BOUT_URL} BOUT-dev && \
    cd BOUT-dev && \
    git checkout ${BOUT_COMMIT} && \
    git submodule update --init --recursive

WORKDIR /home/boutuser/BOUT-dev

# Activate Spack environment and configure BOUT++ with CUDA
USER root
RUN . /spack/share/spack/setup-env.sh && \
    spack env activate /spack-env && \
    spack load cmake && \
    spack load fmt && \
    spack load raja && \
    spack load umpire && \
    spack load netcdf-cxx4 && \
    spack load fftw && \
    chown -R boutuser:boutuser /home/boutuser/BOUT-dev

USER boutuser

# Configure BOUT++ with CUDA support
# Note: cuda_arch=80 corresponds to Ampere architecture (A100)
# Adjust cuda_arch based on your target GPU architecture
RUN . /spack/share/spack/setup-env.sh && \
    spack env activate /spack-env && \
    spack load cmake fmt raja umpire netcdf-cxx4 fftw && \
    cmake -S . -B build \
      -DCMAKE_INSTALL_PREFIX=/opt/bout++ \
      -DCMAKE_BUILD_TYPE=Release \
      -DBOUT_ENABLE_CUDA=ON \
      -DCUDA_ARCH="compute_80,code=sm_80" \
      -DBOUT_ENABLE_RAJA=ON \
      -DBOUT_ENABLE_UMPIRE=ON \
      -DBOUT_USE_SYSTEM_FMT=ON \
      -DBOUT_GENERATE_FIELDOPS=OFF \
      -DBOUT_ENABLE_PYTHON=ON \
      || (cat /home/boutuser/BOUT-dev/build/CMakeFiles/CMake{Output,Error}.log ; exit 1)

# Build BOUT++
RUN . /spack/share/spack/setup-env.sh && \
    spack env activate /spack-env && \
    spack load cmake fmt raja umpire netcdf-cxx4 fftw && \
    cmake --build build -j$(nproc)

# Install BOUT++
USER root
RUN . /spack/share/spack/setup-env.sh && \
    spack env activate /spack-env && \
    spack load cmake && \
    cmake --install build

# Create a helper script to activate the Spack environment
RUN echo '#!/bin/bash\n\
. /spack/share/spack/setup-env.sh\n\
spack env activate /spack-env\n\
spack load cmake fmt raja umpire netcdf-cxx4 fftw\n\
exec "$@"' > /usr/local/bin/bout-env.sh && \
    chmod +x /usr/local/bin/bout-env.sh

USER boutuser
WORKDIR /home/boutuser/BOUT-dev

# Default command activates Spack environment and drops to bash
ENTRYPOINT ["/usr/local/bin/bout-env.sh"]
CMD ["/bin/bash"]
