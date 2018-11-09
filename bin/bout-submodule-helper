#!/bin/sh
# GitHub Releases don't include submodules, and it's not possible to
# use `git submodule update` unless you're in a git repo. So we first
# need to make a git repo, then parse the .gitmodules file and add the
# submodules from that. Finally, we delete the temporary git repo
#
# There is a big danger here, in that the submodules won't be at the
# versions we've committed in the real repo. I can't see a way
set -e

usage()
{
    cat <<EOF
Usage: bout-submodule-helper [DIRECTORY]

Try to download the submodules need to run/test BOUT++

Optional DIRECTORY is the top of the repository
EOF
    exit $1
}

if [ $# -gt 2 ]; then
    usage 1
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --help)
            usage 0
            ;;
        *)
            cd $1
            ;;
    esac
    shift
done

if $(git rev-parse --is-inside-work-tree > /dev/null 2>&1); then
    # Good news! We're in a git repo. We can do the sensible thing
    git submodule update --init --recursive
    exit 0
else
    echo "Not in a git repo: creating temporary repo in order to clone submodules"
    # Better hope README.md exists... If anything goes wrong here, we
    # need to delete the git repo we've made, or it won't be possible
    # to rerun this script.
    git init >/dev/null || (rm -rf ./.git; exit 1)
    git add README.md >/dev/null || (rm -rf ./.git; exit 2)
    git commit -am "DANGER: THIS WAS CREATED AUTOMATICALLY. SUBMODULES MAY BE AT WRONG COMMITS" \
        >/dev/null || (rm -rf ./.git; exit 3)
fi

# The following adapted from https://stackoverflow.com/a/11258810/2043465
git config -f .gitmodules --get-regexp '^submodule\..*\.path$' |
    while read path_key path
    do
        url_key=$(echo $path_key | sed 's/\.path/.url/')
        url=$(git config -f .gitmodules --get "$url_key")
        if [ -d ./$path ]; then
            echo -e "\n\n** WARNING: Directory './$path' already exists! Please delete and try again **\n\n"
            rm -rf ./.git
            exit 4
        fi
        git submodule add $url $path || \
            (echo "Something went wrong downloading the submodules: aborting";
             rm -rf ./.git;
             exit 5)
    done

rm -rf ./.git
