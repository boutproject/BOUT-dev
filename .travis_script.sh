#!/bin/bash
echo "Configuring with $CONFIGURE_OPTIONS"
time ./configure $CONFIGURE_OPTIONS
conf=$?
if test $conf -gt 0
then
    echo
    echo "Printing config.log:"
    echo
    echo
    cat config.log
    echo
    echo "Printing config-build.log:"
    echo
    echo
    cat config-build.log
    exit $conf
fi
export PYTHONPATH=$(pwd)/tools/pylib/:$PYTHONPATH
njobs=4
time make -j ${njobs} || exit
time make -j ${njobs} check-unit-tests || exit
time make -j ${njobs} check-integrated-tests || exit
time make -j ${njobs} check-mms-tests || exit
