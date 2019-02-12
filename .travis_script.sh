#!/bin/bash

#Default flags
COVERAGE=0
UNIT=0
INTEGRATED=0
MMS=0
TESTS=0
MAIN_TARGET=

usage() {
    echo "$0 options are: "
    #Just pull out special comments from this file
    grep "\#\#\#" $0
    exit 1
}

#Handle input flags
while getopts "cuimt:" arg;
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
    t) ### Set target to build
        MAIN_TARGET+=("$OPTARG")
        ;;
	*) ### Show usage message
	    usage
	    ;;
    esac
done

export MAKEFLAGS="-j 2 -k"
echo "Configuring with $CONFIGURE_OPTIONS"
time ./configure $CONFIGURE_OPTIONS MAKEFLAGS="$MAKEFLAGS"
conf=$?
if test $conf -gt 0
then
    RED_FG="\033[031m"
    RESET_FG="\033[039m"
    echo -e $RED_FG
    echo "**************************************************"
    echo "Printing config.log:"
    echo "**************************************************"
    echo -e $RESET_FG
    echo
    cat config.log
    echo
    echo -e $RED_FG
    echo "**************************************************"
    echo "Printing config-build.log:"
    echo "**************************************************"
    echo -e $RESET_FG
    echo
    cat config-build.log
    exit $conf
fi
export PYTHONPATH=$(pwd)/tools/pylib/:$PYTHONPATH

for target in ${MAIN_TARGET[@]}
do
    time make $target
    make_exit=$?
    if [[ $make_exit -gt 0 ]]; then
	make clean > /dev/null
	echo -e $RED_FG
	echo "**************************************************"
	echo "Printing make commands:"
	echo "**************************************************"
	echo -e $RESET_FG
	echo
	make -n $target
	exit $make_exit
    fi
done

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
    time py.test-3 tools/pylib/ || exit
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

    # Upload for codecov. Slightly hacky: for each source file,
    # codecov runs gcov in that directory (i.e. essentially `cd
    # src/sys && gcov utils.cxx`). This is mostly fine, but it seems
    # we also need to run it just in the unit tests directory for each
    # for of the unit test sources (i.e. `cd tests/unit/; gcov
    # sys/test_utils.cxx`). This is the difference between the
    # `-execdir` and `-exec` arguments to `find`
    bash <(curl -s https://codecov.io/bash) -X fix |\
        -a '-r {} +; find ./tests/unit -type f -name '*.gcno' -exec gcov -pbr {} +'

    #For codacy
    bash ./.codacy_coverage.sh
    
fi
