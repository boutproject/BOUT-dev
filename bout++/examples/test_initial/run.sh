#!/bin/bash

make

MPIRUN=mpirun
MD5SUM=md5sum

ntotal=0
npassed=0

############### No ballooning #############

echo "No ballooning"

cd data
rm BOUT.inp
ln -s BOUT.inp_nobal BOUT.inp
cd ..

# Run on 4 processors

echo "   4 processors"
$MPIRUN -np 4 ./test_initial >& log.txt
errmsg=`grep FAILED data/BOUT.log.*` 

if test "$errmsg" = ""; then
    echo "     => TEST PASSED"
    npassed=$[$npassed+1]
else
    echo "     => TEST FAILED"
fi
ntotal=$[$ntotal+1]

idl gather >& idl_output.txt
mv result.nc result_4.nc

# Run on 8 processors
echo "   8 processors"
$MPIRUN -np 8 ./test_initial >& log.txt
errmsg=`grep FAILED data/BOUT.log.*` 

if test "$errmsg" = ""; then
    echo "     => TEST PASSED"
    npassed=$[$npassed+1]
else
    echo "     => TEST FAILED"
fi
ntotal=$[$ntotal+1]

idl gather >& idl_output.txt

# Compare outputs
echo "   Compare results"
hash4=`md5sum result_4.nc | head -c 32`
hash8=`md5sum result.nc | head -c 32`
if test "$hash4" = "$hash8"; then
    echo "     => TEST PASSED"
    npassed=$[$npassed+1]
else
    echo "     => TEST FAILED"
fi
ntotal=$[$ntotal+1]

############### ballooning ##############

echo "Ballooning transform in the core"

cd data
rm BOUT.inp
ln -s BOUT.inp_bal BOUT.inp
cd ..

# Run on 4 processors

echo "   4 processors"
$MPIRUN -np 4 ./test_initial >& log.txt
errmsg=`grep FAILED data/BOUT.log.*` 

if test "$errmsg" = ""; then
    echo "     => TEST PASSED"
    npassed=$[$npassed+1]
else
    echo "     => TEST FAILED"
fi
ntotal=$[$ntotal+1]

idl gather >& idl_output.txt
mv result.nc result_4.nc

# Run on 8 processors
echo "   8 processors"
$MPIRUN -np 8 ./test_initial >& log.txt
errmsg=`grep FAILED data/BOUT.log.*` 

if test "$errmsg" = ""; then
    echo "     => TEST PASSED"
    npassed=$[$npassed+1]
else
    echo "     => TEST FAILED"
fi
ntotal=$[$ntotal+1]

idl gather >& idl_output.txt

# Compare outputs
echo "   Compare results"
hash4=`md5sum result_4.nc | head -c 32`
hash8=`md5sum result.nc | head -c 32`
if test "$hash4" = "$hash8"; then
    echo "     => TEST PASSED"
    npassed=$[$npassed+1]
else
    echo "     => TEST FAILED"
fi
ntotal=$[$ntotal+1]

echo "RESULT: Passed $npassed out of $ntotal tests"
