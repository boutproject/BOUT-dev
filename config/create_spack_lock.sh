#!/usr/bin/env bash

# The following modules are provided by the system
module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cmake/3.14.5
module --ignore-cache load cuda/10.1.243
module --ignore-cache load gcc/8.3.1
module --ignore-cache load sundials/4.1.0
module --ignore-cache load lapack/3.8.0-gcc-4.9.3

if [ ! -d "$SPACK_ROOT/var/spack/environments/BOUT" ] 
then
spack env create BOUT  
fi

spack env activate BOUT

spack add petsc@3.12.3%gcc@8.3.1 +fftw +metis +superlu-dist +hypre +mpi ^hdf5@1.10.1%gcc@8.3.1+cxx+hl+mpi+pic+shared+fortran ^spectrum-mpi@rolling-release%gcc@8.3.1 ^metis%gcc@8.3.1+real64 ^fftw%gcc@8.3.1

spack add netcdf-cxx4@4.3.1%gcc@8.3.1 ^hdf5@1.10.1%gcc@8.3.1+cxx+fortran+hl+mpi+pic+shared ^spectrum-mpi@rolling-release%gcc@8.3.1 

spack add netcdf-c@4.7.3%gcc@8.3.1 ^spectrum-mpi@rolling-release%gcc@8.3.1 

spack add python@3.7.6%gcc@8.3.1
spack add py-setuptools@41.4.0%gcc@8.3.1
spack add py-cftime@1.0.3.4%gcc@8.3.1
spack add py-cython@0.29.14%gcc@8.3.1
spack add py-pybind11@2.5.0%gcc@8.3.1
spack add py-numpy@1.18.2%gcc@8.3.1
spack add py-scipy@1.4.1%gcc@8.3.1

spack add py-netcdf4@1.4.2%gcc@8.3.1  ^spectrum-mpi@rolling-release%gcc@8.3.1 ^hdf5@1.10.1%gcc@8.3.1+cxx+fortran+hl+mpi+pic+shared

spack concretize



