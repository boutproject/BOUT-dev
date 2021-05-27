#!/usr/bin/env bash
cmake --version
cmake . -B build $@ -DCMAKE_INSTALL_PREFIX=$(pwd)/installed
cmake --build build --target build-check -j 2
cd build
ctest --output-on-failure --timeout 300

# Test bout-config basic functionallity
cd ../examples/make-script
PATH=../../build/bin:$PATH make
./test --help
cd -
make install -j 2
cd -
rm test
PATH=../../installed/bin:$PATH make
./test --help
