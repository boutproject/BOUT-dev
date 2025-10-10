#!/usr/bin/env bash

set -ex

cmake --version
cmake -S . -B build $@ -DCMAKE_INSTALL_PREFIX=$(pwd)/installed

if test ".$UNIT_ONLY" = ".YES" ; then
    make -C build build-check-unit-tests
    make -C build check-unit-tests
    exit 0
fi

cmake --build build --target build-check -j 2
cd build
ctest --output-on-failure --timeout 300

export LD_LIBRARY_PATH=/home/runner/local/lib:$LD_LIBRARY_PATH

# Test bout++Config.cmake
cd ../examples/conduction
cmake . -B build -DCMAKE_PREFIX_PATH=../../build
cmake --build build
./build/conduction
cd -

make install -j 2
rm -rf build
# Test installation with CMake
cd ../examples/conduction
cmake . -B build -DCMAKE_PREFIX_PATH=../../installed
cmake --build build
./build/conduction
