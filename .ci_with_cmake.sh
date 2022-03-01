#!/usr/bin/env bash

set -ex

cmake --version
cmake -S . -B build $@ -DCMAKE_INSTALL_PREFIX=$(pwd)/installed

if test ".$UNIT_ONLY" = ".YES" ; then
    make -C build build-check-unit-tests
    make -C build check-unit-tests
else
    cmake --build build --target build-check -j 2
    cd build
    ctest --output-on-failure --timeout 300
fi

export LD_LIBRARY_PATH=/home/runner/local/lib:$LD_LIBRARY_PATH

# Test bout-config basic functionallity
cd ../examples/make-script
PATH=../../build/bin:$PATH bout-config --all
PATH=../../build/bin:$PATH make
./test --help
cd -
make install -j 2
cd -
rm test
PATH=../../installed/bin:$PATH bout-config --all
PATH=../../installed/bin:$PATH make
./test --help
