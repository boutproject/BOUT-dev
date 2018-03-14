#!/bin/bash
echo "Configuring with $CONFIGURE_OPTIONS"
#time ./configure $CONFIGURE_OPTIONS
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

#Default flags
COVERAGE=0
UNIT=0
INTEGRATED=0
MMS=0
TESTS=0

usage() {
    echo "$0 options are: "
    #Just pull out special comments from this file
    grep "\#\#\#" $0}
}
#Handle input flags
while getopts "cuim" arg;
do
    case $arg in
	c) ### Run the coverage-post job tasks
	    COVERAGE=1
	    ;;
	u) ### Run the unit tests
	    UNIT=1
	    TESTS=1
	    ;;
	i) ### Run the integrated tests
	    INTEGRATED=1
	    TESTS=1
	    ;;
	m) ### Run the mms tests
	    MMS=1
	    TESTS=1
	    ;;
	*) ### Show usage message
	    usage
	    ;;
    esac
done

time make $MAIN_TARGET|| exit
if [[ ${TESTS} == 1 ]]
then
    time make build-check || exit
fi

if [[ ${UNIT} == 1 ]]
then
    time make check-unit-tests || exit
fi

if [[ ${INTEGRATED} == 1 ]]
then
    time make check-integrated-tests || exit
fi

if [[ ${MMS} == 1 ]]
then
   time make check-mms-tests || exit
fi

if [[ ${COVERAGE} == 1 ]]
then
   # Ensure that there is a corresponding .gcda file for every .gcno file
   # This is to try and make the coverage report slightly more accurate
   # It still won't include, e.g. any solvers we don't build with though
   find . -name "*.gcno" -exec sh -c 'touch -a "${1%.gcno}.gcda"' _ {} \;
   bash <(curl -s https://codecov.io/bash)
fi
