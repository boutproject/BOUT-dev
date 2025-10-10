#!/bin/bash
set -e

echo "===> Building BOUT-dev CUDA minimal"
cmake -S . -B build \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DBOUT_ENABLE_RAJA=on \
    -DBOUT_ENABLE_UMPIRE=on \
    -DBOUT_ENABLE_CUDA=on \
    -DCMAKE_CUDA_ARCHITECTURES=70 \
    -DCUDA_ARCH=compute_70,code=sm_70 \
    -DBOUT_ENABLE_WARNINGS=off \
    -DBOUT_USE_SYSTEM_FMT=on

pushd build
make -j

echo "===> Building and running blob2d-outerloop"
pushd examples/blob2d-outerloop
make -j
# Check the output using Sim Time and RHS evals. Must be careful splitting the
# regex string in mulitple lines and escaping characters.
if ./blob2d-outerloop | grep -Pzoq "(?s)Sim Time  \|  RHS evals  \| Wall Time \|  Calc    Inv   Comm    I/O   SOLVER\n.*\n"\
"0\.000e\+00          2       .*"\
"5\.000e\+01         53       .*"\
"1\.000e\+02         17       .*"\
"1\.500e\+02         27       .*"; then
    echo "Sim Time and RHS evals match"
else
    echo "Sim Time and RHS evals DO NOT match"
    exit 1
fi
popd

echo "===> Building and running elm-pb-outerloop"
pushd examples/elm-pb-outerloop
make -j
if ./elm_pb_outerloop | grep -Pzoq "(?s)Sim Time  \|  RHS evals  \| Wall Time \|  Calc    Inv   Comm    I/O   SOLVER\n.*\n"\
"0\.000e\+00          2       .*"\
"1\.000e\+00         44       .*"\
"2\.000e\+00         37       .*"\
"3\.000e\+00         37       .*"\
"4\.000e\+00         37       .*"\
"5\.000e\+00         30       .*"\
"6\.000e\+00         31       .*"\
"7\.000e\+00         31       .*"\
"8\.000e\+00         25       .*"\
"9\.000e\+00         21       .*"\
"1\.000e\+01         24       .*"\
"1\.100e\+01         19       .*"\
"1\.200e\+01         25       .*"\
"1\.300e\+01         25       .*"\
"1\.400e\+01         25       .*"\
"1\.500e\+01         25       .*"\
"1\.600e\+01         25       .*"\
"1\.700e\+01         25       .*"\
"1\.800e\+01         25       .*"\
"1\.900e\+01         20       .*"\
"2\.000e\+01         29       .*"\
"2\.100e\+01         29       .*"\
"2\.200e\+01         29       .*"\
"2\.300e\+01         29       .*"\
"2\.400e\+01         29       .*"\
"2\.500e\+01         29       .*"\
"2\.600e\+01         29       .*"\
"2\.700e\+01         22       .*"\
"2\.800e\+01         29       .*"\
"2\.900e\+01         29       .*"\
"3\.000e\+01         29       .*"\
"3\.100e\+01         29       .*"\
"3\.200e\+01         29       .*"\
"3\.300e\+01         32       .*"\
"3\.400e\+01         25       .*"\
"3\.500e\+01         33       .*"\
"3\.600e\+01         33       .*"\
"3\.700e\+01         39       .*"\
"3\.800e\+01         31       .*"\
"3\.900e\+01         31       .*"\
"4\.000e\+01         36       .*"; then
    echo "Sim Time and RHS evals match"
else
    echo "Sim Time and RHS evals DO NOT match"
    exit 1
fi
popd

popd