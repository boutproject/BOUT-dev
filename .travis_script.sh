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
export MAKEFLAGS="-j 2 -k"
time make $MAIN_TARGET|| exit
time make build-check || exit
time make check-unit-tests || exit
time make check-integrated-tests || exit
time make check-mms-tests || exit
