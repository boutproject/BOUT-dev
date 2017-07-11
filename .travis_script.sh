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
time make || exit
time make check || exit
export PYTHONPATH=$(pwd)/tools/pylib/:$PYTHONPATH
cd ./examples
time ./test_suite_make && ./test_suite
