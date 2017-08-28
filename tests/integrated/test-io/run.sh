#!/bin/bash

make

./test_io

benchmark=`md5sum data/benchmark.out.0.nc | head -c 32`
output=`md5sum data/test_io.out.0.nc | head -c 32`

echo "== Checksums =="
echo $benchmark
echo $output
echo ""

if test "$benchmark" = "$output"; then
    echo "=> TEST PASSED"
else
    echo "=> TEST FAILED"
fi

