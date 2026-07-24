#!/usr/bin/env bash

set -ex

cmake --version
cmake -S . -B build $@ -DCMAKE_INSTALL_PREFIX=$(pwd)/installed

cmake --build build --target help
cmake --build build --target build-check -j 2
cd build/tests/MMS/spatial/fci/
PYTHONPATH=../../../..//tools/pylib/: ./runtest
