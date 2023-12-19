#!/bin/bash

set -e

if test $BUILD_ADIOS2 ; then
    if [[ ! -d $HOME/local/adios/include/adios2.h ]] || test $1 ; then
	echo "****************************************"
	echo "Building ADIOS2"
	echo "****************************************"

	branch=${1:-release_29}
    if [ ! -d adios2 ]; then
    	git clone -b $branch https://github.com/ornladios/ADIOS2.git adios2 --depth=1
    fi

	pushd adios2
    rm -rf build
    mkdir -p build
    pushd build

    cmake .. \
        -DCMAKE_INSTALL_PREFIX=$HOME/local \
        -DADIOS2_USE_MPI=ON \
        -DADIOS2_USE_Fortran=OFF \
        -DADIOS2_USE_Python=OFF \
        -DADIOS2_BUILD_EXAMPLES=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBUILD_TESTING=OFF \
        -DADIOS2_USE_SST=OFF \
        -DADIOS2_USE_MGARD=OFF \
        -DADIOS2_USE_HDF5=OFF \
        -DADIOS2_USE_BZip2=OFF \
        -DADIOS2_USE_Blosc2=OFF \
        -DADIOS2_USE_SZ=OFF \
        -DADIOS2_USE_ZFP=OFF \
        -DADIOS2_USE_DAOS=OFF \
        -DADIOS2_USE_UCX=OFF \
        -DADIOS2_USE_LIBPRESSIO=OFF \
        -DADIOS2_USE_Sodium=OFF \
        -DADIOS2_USE_ZeroMQ=OFF \
        -DADIOS2_USE_MHS=OFF \
        -DADIOS2_USE_DataMan=OFF

	make -j 4 && make install
	popd

	echo "****************************************"
	echo " Finished building ADIOS2"
	echo "****************************************"

    else

	echo "****************************************"
	echo " ADIOS2 already installed"
	echo "****************************************"
    fi
else
    echo "****************************************"
    echo " ADIOS2 not requested"
    echo "****************************************"
fi
