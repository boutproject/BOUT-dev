#!/bin/bash

set -e

if [[ ! -d $HOME/local/include/sundials ]]; then
    echo "****************************************"
    echo "Building SUNDIALS"
    echo "****************************************"
    sundials_ver=4.1.0
    wget https://computation.llnl.gov/projects/sundials/download/sundials-${sundials_ver}.tar.gz
    tar xvf sundials-${sundials_ver}.tar.gz
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
