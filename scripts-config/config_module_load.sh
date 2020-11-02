#!/usr/bin/env bash
# The following modules are provided by the system
module --ignore-cache load spectrum-mpi/rolling-release
module --ignore-cache load cmake/3.14.5
module --ignore-cache load cuda/10.1.243
module --ignore-cache load gcc/8.3.1
module --ignore-cache load sundials/4.1.0
#module --ignore-cache load sundials-5.1.0-gcc-8.3.1-g7h45nh
module --ignore-cache  lapack/3.9.0-gcc-7.3.1

# The following modules are installed by you using spack
module --ignore-cache load openblas-0.3.7-gcc-8.3.1-3izwwk6
module --ignore-cache load fftw-3.3.8-gcc-8.3.1-vlusxnt
module --ignore-cache load  hdf5-1.10.1-gcc-8.3.1-xkc527f
module --ignore-cache load metis-5.1.0-gcc-8.3.1-wbsy3pr
module --ignore-cache load parmetis-4.0.3-gcc-8.3.1-7lwyokt
module --ignore-cache load superlu-dist-6.1.1-gcc-8.3.1-scsdvi6
module --ignore-cache load  python-3.7.6-gcc-8.3.1-z4dxpor
module --ignore-cache load netcdf-4.7.2-gcc-8.3.1-nb5sckv # auto installed as part of netcdf-cxx4 install
module --ignore-cache load zlib-1.2.11-gcc-8.3.1-drkbhfs
module --ignore-cache load netcdf-cxx4-4.3.1-gcc-8.3.1-x2rnmfd
module --ignore-cache load petsc-3.12.3-gcc-8.3.1-uasdxmu 

module --ignore-cache load py-cftime-1.0.3.4-gcc-8.3.1-qfrdsh5
module --ignore-cache load py-cython-0.29.14-gcc-8.3.1-axnyudk
module --ignore-cache load py-pybind11-2.5.0-gcc-8.3.1-66yacp2
module --ignore-cache load py-numpy-1.17.3-gcc-8.3.1-blxr3cn
module --ignore-cache load py-scipy-1.4.1-gcc-8.3.1-m2jgedc
module --ignore-cache load py-netcdf4-1.4.2-gcc-8.3.1-xtg4tg2

