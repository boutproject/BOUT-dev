#!/bin/dash
cache=.findall.cache

if test -f $cache && ! test $(find $cache -mtime +1)
then
    cat $cache
    exit 0
fi

get() {
    root=$(pwd)
    cd $1 || exit 1
    makefile=$(perl -p -e 's/\\\n//' makefile)
    for f in $(echo "$makefile" | grep ^[[:space:]]*SOURCEC[[:space:]]*= |cut -d= -f 2)
    do
        echo $1/$f
    done

    for d in $(echo "$makefile" | grep ^[[:space:]]*DIRS[[:space:]]*= |cut -d= -f 2)
    do
        cd $root
        get $1/$d
    done
}

(get .)> $cache

cat $cache
