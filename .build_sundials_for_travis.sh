#!/bin/bash

set -e

if [[ ! -d $HOME/local/include/sundials ]]; then
    echo "****************************************"
    echo "Building SUNDIALS"
    echo "****************************************"
    sundials_ver=4.1.0
    tarball_name=sundials-${sundials_ver}.tar.gz
    llnl_url=https://computation.llnl.gov/projects/sundials/download/${tarball_name}
    github_url=https://github.com/LLNL/sundials/releases/download/v${sundials_ver}/${tarball_name}
    wget ${llnl_url} \
        || wget ${github_url} \
        || (echo "Could not download SUNDIALS! Aborting"; exit 1)
    tar xvf ${tarball_name}
    mkdir -p sundials-${sundials_ver}/build && cd sundials-${sundials_ver}/build
    cmake -DCMAKE_INSTALL_PREFIX="$HOME/local" \
          -DEXAMPLES_INSTALL=off \
          -DMPI_ENABLE=on \
          -DOPENMP_ENABLE=off \
          -DBUILD_CVODES=off \
          -DBUILD_IDAS=off \
          -DBUILD_KINSOL=off \
          -DBUILD_TESTING=off \
          -DMPI_C_COMPILER="$(command -v mpicc)" \
          -DMPI_CXX_COMPILER="$(command -v mpic++)" \
          -DMPIEXEC_EXECUTABLE="$(command -v mpiexec)" \
          ..
    make && make install
    cd "${TRAVIS_BUILD_DIR}"
    echo "****************************************"
    echo "Finished building SUNDIALS"
    echo "****************************************"
else
    echo "****************************************"
    echo "SUNDIALS already installed"
    echo "****************************************"
fi
