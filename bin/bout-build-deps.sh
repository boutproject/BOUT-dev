#!/usr/bin/env bash

PREFIX=${PREFIX:-$HOME/.local/bout-deps}
BUILD=${BUILD:-$HOME/soft/bout-deps}
MKFLAGS=${MKFLAGS:-"-j 16"}
BOUT_TOP=${BOUT_TOP:-$(pwd)}

HDF5VER=${HDF5VER:-1.12.0}
NCVER=${NCVER:-4.7.4}
NCCXXVER=${NCCXXVER:-4.3.1}
FFTWVER=${FFTWVER:-3.3.9}
SUNVER=${SUNVER:-5.7.0}
PETSCVER=${PETSCVER:-3.15.0}


HDF5FLAGS=${HDF5FLAGS:-}
NCFLAGS=${NCFLAGS:-}
NCCXXFLAGS=${NCCXXFLAGS:-}
FFTWFLAGS=${FFTWFLAGS:---enable-avx512 --enable-avx-128-fma}
SUNFLAGS=${SUNFLAGS:-}
PETSCFLAGS=${PETSCFLAGS:-}

CHECK=${CHECK:-yes} # Run checks: 'yes' to run - 'no' to skip
HDF5CHECK=${HDF5CHECK:-${CHECK}}
NCCHECK=${NCCHECK:-${CHECK}}
NCCXXCHECK=${NCCXXCHECK:-${CHECK}}
FFTWCHECK=${FFTWCHECK:-${CHECK}}
SUNCHECK=${SUNCHECK:-${CHECK}}
PETSCCHECK=${PETSCCHECK:-${CHECK}}


help() {
echo The following options are available:
echo
grep '=${.*:-' $0 | grep -v grep| sed -E "s/.*\\\$\{([A-Z0-9a-z_]*):-(.*)\}/Option \$\1 with default '\2'/"
cat <<EOF

Set any of the options by setting the respective environmente
variables, e.g.

   FFTWVER=3.3.8 $0

if you want to set FFTW to an older version.

For more options edit the script or use as a starting point for your
own scripts.  Especially feel free to comment some sections towards
the end of the file, to disable building some parts, search for
"EDIT BELOW" to where to start.

EOF

# Exit of 0 if -h or --help - else 1
test ".$1" = .-h || test ".$1" = .--help
exit
}

setup() {
    mkdir -p $PREFIX
    mkdir -p $PREFIX/lib
    mkdir -p $PREFIX/bin
    if test -e $PREFIX/lib64 ; then
        if ! test -L $PREFIX/lib64 ; then
            echo "$PREFIX/lib64 exists and is not a symlink to lib - aborting"
            exit 1
        fi
    else
        ln -s lib $PREFIX/lib64
    fi
    mkdir -p $BUILD
}

hdf5() {
    cd $BUILD
    wget -c https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-${HDF5VER}/src/hdf5-${HDF5VER}.tar.bz2 || :
    tar -xvf hdf5-$HDF5VER.tar.bz2
    cd hdf5-${HDF5VER}
    ./configure --prefix $PREFIX --enable-build-mode=production $HDF5FLAGS
    make $MKFLAGS
    if test $HDF5CHECK = yes ; then
	make check
    fi
    make install
}

netcdf() {
    cd $BUILD
    wget -c https://github.com/Unidata/netcdf-c/archive/v$NCVER/netcdf-$NCVER.tar.gz || :
    tar -xf netcdf-$NCVER.tar.gz
    cd netcdf-c-$NCVER
    CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib/" ./configure --prefix=$PREFIX $NCFLAGS
    make $MKFLAGS
    if test $NCCHECK = yes ; then
	make check
    fi
    make install
}

nccxx() {
    cd $BUILD
    wget -c ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx4-$NCCXXVER.tar.gz || :
    tar -xf netcdf-cxx4-$NCCXXVER.tar.gz
    cd netcdf-cxx4-$NCCXXVER
    CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib/" ./configure --prefix=$PREFIX $NCCXXFLAGS
    make $MKFLAGS
    if test $NCCXXCHECK = yes ; then
	make check
    fi
    make install
}

fftw() {
    cd $BUILD
    wget -c http://www.fftw.org/fftw-$FFTWVER.tar.gz || :
    tar -xf fftw-$FFTWVER.tar.gz
    cd fftw-$FFTWVER
    ./configure --prefix $PREFIX --enable-shared --enable-sse2 --enable-avx --enable-avx2  $FFTWFLAGS
    make $MKFLAGS
    if test $FFTWCHECK = yes ; then
	make check
    fi
    make install
}

sundials() {
    cd $BUILD
    wget -c https://github.com/LLNL/sundials/archive/v$SUNVER/sundials-$SUNVER.tar.gz || :
    tar -xvf sundials-$SUNVER.tar.gz
    cd sundials-$SUNVER
    mkdir -p build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DMPI_ENABLE=ON .. $SUNFLAGS
    make $MKFLAGS
    if test $SUNCHECK = yes ; then
	make test
    fi
    make install
}

error () {
    echo "$@" >&2
    echo "$@"
    exit 1
}

petsc() {
    test -z $PETSC_DIR || error "\$PETSC_DIR is set ($PETSC_DIR) - please unset"
    test -z $PETSC_ARCH || error "\$PETSC_ARCH is set ($PETSC_ARCH) - please unset"
    cd $BUILD
    wget -c https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$PETSCVER.tar.gz || :
    tar -xf petsc-$PETSCVER.tar.gz
    cd petsc-$PETSCVER
    unset PETSC_DIR
    ./configure COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3" --with-batch --known-mpi-shared-libraries=1 --with-mpi-dir=$OPENMPI_HOME --download-fblaslapack \
        --known-64-bit-blas-indices=0 --download-hypre --with-debugging=0 --prefix=$PREFIX $PETSCFLAGS
    make $MKFLAGS
    if test $PETSCCHECK = yes ; then
	make check
    fi
    make install
}

submod() {
    cd $BOUT_TOP
    git submodule update --init --recursive
}


info() {
    set +x
    cat <<EOF
As an alternative to using environment module you can set the following variables directly:

  export BOUT_DEP=$PREFIX
  export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH

Run configure with:

  ./configure --with-netcdf=\$BOUT_DEP --with-sundials=\$BOUT_DEP --with-fftw=\$BOUT_DEP --with-petsc=\$BOUT_DEP

or

  ./configure --with-netcdf=$BOUT_DEP --with-sundials=$BOUT_DEP --with-fftw=$BOUT_DEP --with-petsc=$BOUT_DEP

EOF
}

envmodule() {
    set +x
    found=$(echo "$MODULEPATH": | while read -d : d ;
    do
        if test -w $d ; then
            echo $d
            break
        fi
    done)
    reload=
    if test ".$found" = . ; then
        MODDIR=${MODDIR:-$HOME/.module} # Path to add to $MODULEPATH
        msg=
        if test -w $MODDIR && test -d $MODDIR ; then
            # All good
            found=$MODDIR
        elif test -e $MODDIR ; then
            error $MODDIR exists but is not a writeable dir - please set \$MODDIR to something else
        else
            found=$MODDIR
            msg=" create the folder and"
        fi
        echo -n "Do you want to$msg add $MODDIR to \$MODULEPATH to your .bashrc? [Y/n]"
        read ans
        test -z "$ans" || test "$ans" = y || test "$ans" = Y || return
        mkdir -p $MODDIR
        cat <<EOF >> ~/.bashrc
# Add local module path, added by $0
MODULEPATH=\$MODULEPATH:$MODDIR
EOF
        found="$MODDIR"
        reload=yes
    fi
    ext0=bout-dep/$(date +"%Y-%m-%d")
    ext=$ext0
    of="$found/$ext"
    if test -e $of ; then
        of0=$of
        i=0
        while test -e $of~$i ; do i=$((i+1)) ; done
        of="$of~$i"
        ext="$ext~$i"
    fi
    test -e $MODDIR/bout-dep || mkdir $MODDIR/bout-dep
    echo "#%Module 1.0
#
#  BOUT++ dependency module for use with 'environment-modules' package
#


# Only allow one bout-dep module to be loaded at a time
conflict bout-dep
# Require all modules that where loaded at generation time
prereq $(echo $LOADEDMODULES | tr : \ )


setenv        BOUT_DEP         $PREFIX
prepend-path  PATH             $PREFIX/bin
prepend-path  LD_LIBRARY_PATH  $PREFIX/lib
# Avoid potential conflicts
unsetenv      PETSC_DIR
unsetenv      PETSC_ARCH
""" > $of
    test $of0 && diff $of0 $of -q && rm $of && of=$dup && ext=$ext0
    test $reload && msg="
Re-login in or source ~/.bashrc
" || msg=
    cat <<EOF
####################################################################
###                  Summary and Info follows!                   ###
####################################################################
$msg
Activate the dependency module by running

  module load $ext

EOF
}

test ".$1" = . || help $@

# - EDIT BELOW - EDITBELOW - editbelow - "EDIT BELOW" -
# Feel free to edit below

# Comment this if you want to ignore errors
set -e
# Comment this if you want to have less output
set -x

## Setup folders and links
setup
## Build and install hdf5
hdf5
## Build and install netcdf
netcdf
## Build and install C++ interface for netcdf
nccxx
## Build and install FFTW
fftw
## Build and install Sundials
sundials
## Build and install PETSc
petsc
## Download BOUT++ submodules
submod
## Create a moduleinfo
envmodule
## Print infos on how to proceed (assumes all other steps have been run)
info
