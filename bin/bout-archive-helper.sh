#!/bin/bash
#
# A simple wrapper around git_archive_all.sh

set -e

readonly PROGRAM=$(basename $0)

function usage () {
    echo "usage: $PROGRAM [--format FORMAT] [--quiet] tag|commit"
    echo "       $PROGRAM --help|-h"
    echo
    echo "    Helper utility for running git_archive_all.sh"
    echo
    echo "    If no --format is given, use tar.gz"
    echo
    echo "    If tag or commit is not valid, uses HEAD instead, although the tarball"
    echo "    is still called 'BOUT++-<tag>.tar.gz"
}

git_archive_all_bin=externalpackages/git-archive-all.sh/git-archive-all.sh

# Defaults
format=tar.gz
verbose="--verbose"

readonly NOT_ENOUGH_ARGS=2

if [[ $# == 0 ]]; then
    usage
    exit $NOT_ENOUGH_ARGS
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        --format )
            shift
            format="$1"
            shift
            ;;
        --quiet | -q )
            verbose=""
            shift
            ;;
        --help | -h )
            usage
            exit
            ;;
        * )
            break
            ;;
    esac
done

# Tag or commit to attempt to bundle
treeish=$1

# Tarball version name
bout_name=$treeish

# Common prefix for the tarball file name and internal prefix
bout_prefix="BOUT++-${bout_name}"

if ! $(git rev-parse $treeish >/dev/null 2>&1); then
    if [[ -n $verbose ]]; then
        echo "Warning! $treeish is not a valid tag; using HEAD instead"
    fi
    treeish="HEAD"
fi

if [[ -n $verbose ]]; then
    echo "creating tarball for $bout_prefix using $treeish"
fi

$git_archive_all_bin ${verbose} \
                     --format ${format} \
                     --tree-ish $treeish \
                     --prefix ${bout_prefix}/ \
                      ${bout_prefix}.${format}
