#!/bin/bash

set -ex

#Default flags
COVERAGE=0
UNIT=0
INTEGRATED=0
MMS=0
TESTS=0
MAIN_TARGET=
UPDATE_SCRIPT=0
CONFIGURE_SHELL=

usage() {
    echo "$0 options are: "
    #Just pull out special comments from this file
    grep "\#\#\#" $0
    exit 1
}

#Handle input flags
while getopts "cuimt:5s:" arg;
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
            MAIN_TARGET="$OPTARG"
            ;;
        5) ### Run the update to version 5 script
            UPDATE_SCRIPT=1
            ;;
	s) ### Use specific shell to configure
	    CONFIGURE_SHELL="$OPTARG"
	    ;;
        *) ### Show usage message
	    usage
	    ;;
    esac
done

./.build_sundials_for_travis.sh

if test $UPDATE_SCRIPT -gt 0
then
    # Make sure the header list is up to date
    if ! diff bin/bout_4to5_header_file_list <(cd include/;ls *xx|grep -v ^bout.hxx|sort)
    then
	echo "Some header files changed."
	echo "Please update the list by running:"
	echo "(cd include/;ls *xx|grep -v ^bout.hxx|sort) > bin/bout_4to5_header_file_list"
	echo "And commit the updated file."
	exit 1
    fi

    bin/bout_4to5 -f
fi

export MAKEFLAGS="-j 2 -k"
echo "****************************************"
echo "Configuring with $CONFIGURE_OPTIONS"
echo "****************************************"
conf=0
time $CONFIGURE_SHELL ./configure $CONFIGURE_OPTIONS MAKEFLAGS="$MAKEFLAGS" || conf=$?
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
    make_exit=0
    time make $target || make_exit=$?
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
    time make build-check
fi

if [[ ${UNIT} == 1 ]]
then
    time make check-unit-tests
fi

if [[ ${INTEGRATED} == 1 ]]
then
    time make check-integrated-tests
    time py.test-3 tools/pylib/
fi

if [[ ${MMS} == 1 ]]
then
   time make check-mms-tests
fi

if [[ ${COVERAGE} == 1 ]]
then
    # Ensure that there is a corresponding .gcda file for every .gcno file
    # This is to try and make the coverage report slightly more accurate
    # It still won't include, e.g. any solvers we don't build with though
    find . -name "*.gcno" -exec sh -c 'touch -a "${1%.gcno}.gcda"' _ {} \;

    # Use lcov to generate a report, upload it to codecov.io
    make code-coverage-capture
    bash <(curl -s https://codecov.io/bash) -f bout-coverage.info

    #For codacy
    bash ./.codacy_coverage.sh
fi
