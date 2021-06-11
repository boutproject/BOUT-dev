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
the end of the file, to disable building some parts.

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
    make install
}

netcdf() {
    cd $BUILD
    wget -c https://github.com/Unidata/netcdf-c/archive/v$NCVER/netcdf-$NCVER.tar.gz || :
    tar -xf netcdf-$NCVER.tar.gz
    cd netcdf-c-$NCVER
    CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib/" ./configure --prefix=$PREFIX $NCFLAGS
    make $MKFLAGS
    make install
}

nccxx() {
    cd $BUILD
    wget -c ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx4-$NCCXXVER.tar.gz || :
    tar -xf netcdf-cxx4-$NCCXXVER.tar.gz
    cd netcdf-cxx4-$NCCXXVER
    CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib/" ./configure --prefix=$PREFIX $NCCXXFLAGS
    make $MKFLAGS
    make install
}

fftw() {
    cd $BUILD
    wget -c http://www.fftw.org/fftw-$FFTWVER.tar.gz || :
    tar -xf fftw-$FFTWVER.tar.gz
    cd fftw-$FFTWVER
    ./configure --prefix $PREFIX --enable-shared --enable-sse2 --enable-avx --enable-avx2  $FFTWFLAGS
    make $MKFLAGS
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
    make install
}

submod() {
    cd $BOUT_TOP
    git submodule update --init --recursive
}


info() {
set +x
    echo "Put this in a file in your module path"
    echo "#---------------------------"
    echo "#%Module 1.0
#
#  BOUT++ module for use with 'environment-modules' package


# Only allow one bout-dep module to be loaded at a time
conflict bout-dep
# Require all modules that where loaded at generation time
prereq $(echo $LOADEDMODULES | tr : \ )


setenv        BOUT_DEP         $PREFIX
prepend-path  PATH             $PREFIX/bin
prepend-path  LD_LIBRARY_PATH  $PREFIX/lib
"
    echo "#---------------------------"
    echo Run configure with:
    echo ./configure --with-netcdf=\$BOUT_DEP --with-sundials=\$BOUT_DEP --with-fftw=\$BOUT_DEP --with-petsc=\$BOUT_DEP
}

# Uncomment this if want to use it as a script to detect errors
set -ex
test ".$1" = . || help $@

setup
hdf5
netcdf
nccxx
fftw
sundials
petsc
submod
info
