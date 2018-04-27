#!/bin/sh

cache=.findall.cache

#rm $cache

if test -f $cache && ! test $(find $cache -mtime +1)
then
    cat $cache
    exit 0
fi

get() {
    ## Version 2:
    ## use make to generate a recursive list
    FINDALLCXX='yes' CUDIR=. $MAKE -r -R listfiles -s --no-print-directory
    ## Version 1:
    ## use make to parse the makefiles, and grep for the list of files and directories
    # all=$(  -r  -p  -q -C $1|grep -e '^SOURCEC\|^CURDIR\|^DIRS')
    # for f in $(echo "$all"|grep SOURCE|cut -f2 -d=)
    # do
    #     echo $1/$f
    # done
    # #echo "$all"
    # DIRS=$(echo "$all"|grep CURDIR -B 1| tail -n2 |grep DIRS|cut -f2 -d=)
    # for d in $DIRS
    # do
    #     get $1/$d
    # done
    ## Version 0:
    ## Parse the makefile with perl + grep and find DIRS and SOURCEC
    # makefile=$(perl -p -e 's/\\\n//' makefile)
    # for f in $(echo "$makefile" | grep ^[[:space:]]*SOURCEC[[:space:]]*= |cut -d= -f 2)
    # do
    #     echo $1/$f
    # done
    # for d in $(echo "$makefile" | grep ^[[:space:]]*DIRS[[:space:]]*= |cut -d= -f 2)
    # do
    #     cd $root
    #     get $1/$d
    # done
}

(get .)> $cache

cat $cache
