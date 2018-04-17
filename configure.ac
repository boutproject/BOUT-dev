# Copyright 2010 B D Dudson, S Farley
#
# Contact Ben Dudson, bd512@york.ac.uk
#
# This file is part of BOUT++.
#
# BOUT++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BOUT++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################
#
# Process this file with autoreconf to produce a configure script.
#
#     $ autoreconf -if 
#
# Changelog:
#
# 2010-03-09 Ben Dudson <bd512@york.ac.uk>
#    * Changing to always require FFTW (removing NR routines)
# 2015-08-08 David Schwörer <schword2@mail.dcu.ie>
#    * Searching for libs in lib and lib64
#

AC_PREREQ(2.59)
AC_INIT(BOUT++, 4.1.2, bd512@york.ac.uk)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH(netcdf,       [AS_HELP_STRING([--with-netcdf],
        [Enable support for netCDF files])],,[])
AC_ARG_WITH(pnetcdf,      [AS_HELP_STRING([--with-pnetcdf],
        [Set path to Parallel NetCDF library])],,[])
AC_ARG_WITH(ida,          [AS_HELP_STRING([--with-ida=/path/to/ida],
        [Use the SUNDIALS IDA solver])],,[])
AC_ARG_WITH(cvode,        [AS_HELP_STRING([--with-cvode],
        [Use the SUNDIALS CVODE solver])],,[])
AC_ARG_WITH(sundials,     [AS_HELP_STRING([--with-sundials],
        [Use CVODE and IDA])],,[])
AC_ARG_WITH(fftw,         [AS_HELP_STRING([--with-fftw],
        [Set directory of FFTW3 library])],,[])
AC_ARG_WITH(lapack,       [AS_HELP_STRING([--with-lapack],
        [Use the LAPACK library])],,[])
AC_ARG_WITH(petsc,        [AS_HELP_STRING([--with-petsc],
        [Enable PETSc interface])],,[])
AC_ARG_WITH(slepc,        [AS_HELP_STRING([--with-slepc],
        [Enable SLEPc interface])],,[])
AC_ARG_WITH(pvode,        [AS_HELP_STRING([--with-pvode],
        [Build and enable PVODE 98 (DEFAULT)])],,[])
AC_ARG_WITH(mumps,        [AS_HELP_STRING([--with-mumps],
        [Link with MUMPS library for direct matrix inversions])],,[])
AC_ARG_WITH(arkode,       [AS_HELP_STRING([--with-arkode],
        [Use the SUNDIALS ARKODE solver])],,[])
AC_ARG_WITH(scorep,       [AS_HELP_STRING([--with-scorep],
        [Enable support for scorep based instrumentation])],,[with_scorep=no])

dnl --with-hdf5 flags are set in AX_LIB_{PARALLEL}HDF5

AC_ARG_ENABLE(checks,       [AS_HELP_STRING([--enable-checks=no/1/2/3],
        [Set run-time checking level])],,[])
AC_ARG_ENABLE(signal,       [AS_HELP_STRING([--disable-signal],
        [Disable SEGFAULT handling])],,[])
AC_ARG_ENABLE(color,        [AS_HELP_STRING([--disable-color],
        [Disable -c option to color output])],,[])
AC_ARG_ENABLE(track,        [AS_HELP_STRING([--enable-track],
        [Enable variable tracking])],,[])
AC_ARG_ENABLE(debug,        [AS_HELP_STRING([--enable-debug],
        [Enable all debugging flags])],,[])
AC_ARG_ENABLE(optimize,     [AS_HELP_STRING([--enable-optimize=1/2/3],
        [Enable optimization])],,[])
AC_ARG_ENABLE(sigfpe,       [AS_HELP_STRING([--enable-sigfpe],
        [Enable FloatingPointExceptions])],,[])
AC_ARG_ENABLE(backtrace,    [AS_HELP_STRING([--disable-backtrace],
        [Disable function backtrace])],,[enable_backtrace=maybe])
AC_ARG_ENABLE(shared,       [AS_HELP_STRING([--enable-shared],
        [Enable building bout++ into an shared object])],,[enable_shared=no])
AC_ARG_ENABLE(openmp,       [AS_HELP_STRING([--enable-openmp],
        [Enable building with OpenMP support])],,[enable_openmp=no])
AC_ARG_ENABLE(pvode_openmp,       [AS_HELP_STRING([--enable-pvode-openmp],
        [Enable building PVODE with OpenMP support])],,[enable_pvode_openmp=no])

AC_ARG_VAR(EXTRA_INCS,[Extra compile flags])
AC_ARG_VAR(EXTRA_LIBS,[Extra linking flags])

file_formats=""  # Record which file formats are being supported

# Delete the build log from last time
rm config-build.log

AC_ARG_VAR(CXXFLAGS,[Extra compile flags])
AC_ARG_VAR(LDFLAGS,[Extra linking flags])
AC_ARG_VAR(LDLIBS,[Extra linking libraries])

AC_SUBST(MKDIR_P)
AC_SUBST(EXTRA_INCS)
AC_SUBST(EXTRA_LIBS)

# Adding variables for additional sources
AC_SUBST(PRECON_SOURCE)

# We're using C++
AC_LANG(C++)

#############################################################
# Checks for programs
#############################################################

# Search for MPI compiler; fail if not found
AX_PROG_CXX_MPI([], [], [
  AC_MSG_ERROR([*** An MPI compiler is required. You might need to set MPICXX correctly.])
])

# Utility programs
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_INSTALL

# Set MAKE to gmake if possible, otherwise make
AC_CHECK_PROG(MAKE, gmake, gmake, make)

AC_PROG_RANLIB

AC_SUBST(ARFLAGS)

ARFLAGS=''
for flag in cruU cru
do
  echo 1 > artest1
  ar $flag artest artest1 &&ar $flag artest artest1
  arexit=$?
  rm -f artest1 artest
  if test $arexit -eq 0
  then
    ARFLAGS="$flag"
    break;
  fi
done
test -z $ARFLAGS && AC_MSG_ERROR([Failed to find suitable flags for ar])

# Check for and enable C++11 support
# Error if not supported
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])

#############################################################
# STD Library functions
#############################################################

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h strings.h])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_VPRINTF

# Check for OpenMP support
: ${enable_openmp=no}  # Disable by default
AC_OPENMP
HAS_OPENMP=$enable_openmp

#############################################################
# Code coverage using gcov
#
# Mutally exclusive with optimisation, therefore needs to come first
# so we can turn off optimisation if coverage is enabled
#############################################################

AX_CODE_COVERAGE()
AS_IF([test "x$enable_code_coverage" = "xyes"],
[
    AS_IF([test "x$enable_optimize"], [
        AC_MSG_WARN([Code coverage clashes with optimisations, disabling optimisations])
        enable_optimize="0"
    ])
    COVERAGE_FLAGS="--coverage --no-inline"
    LDFLAGS="$LDFLAGS --coverage"
], [
    COVERAGE_FLAGS=
])
AC_SUBST([COVERAGE_FLAGS])

#############################################################
# General Options
#############################################################

OPT_FLAGS=""
if test "$enable_debug" != ""
then
  echo "Enabling all debug options"
  enable_checks="3"
  # use -Og with available, otherwise fall back to -O0
  OPT_FLAGS="-g -O0 -Og -fno-inline -hipa1"
elif test $enable_optimize ; then
  optimize=$enable_optimize
  test $optimize -eq 0 2>/dev/null
  if test $? -eq 2 ; then
    if test $optimize = "default"
    then
      OPT_FLAGS="-O"
    elif test $optimize = "fast"
    then
      OPT_FLAGS="-Ofast -fno-finite-math-only -march=native -funroll-loops"
      DISABLE_CHECK=probably
    else
      AC_MSG_ERROR([unrecognized option: --enable-optimize=$optimize])
    fi
  else
    if test $optimize -gt 3
    then
      OPT_FLAGS="-Ofast -fno-finite-math-only -march=native -funroll-loops"
      DISABLE_CHECK=probably
    elif test $optimize -gt 2
    then
      OPT_FLAGS="-O3 -march=native -funroll-loops"
      DISABLE_CHECK=probably
    elif test $optimize -gt 1
    then
      OPT_FLAGS="-O2 -march=native"
    elif test $optimize -ge 0
    then
      OPT_FLAGS="-O$optimize"
    else
      AC_MSG_ERROR([unrecognized option: --enable-optimize=$optimize])
    fi
  fi
else
  OPT_FLAGS="-O"
fi
for flag in $OPT_FLAGS
do
  CXXFLAGS_save=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $flag"
  AC_MSG_CHECKING([for flag $flag])
  AC_TRY_COMPILE([],[],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no);CXXFLAGS=$CXXFLAGS_save])
done

# disable checks if optimization > 2 is used
if test -z $enable_checks && test .$DISABLE_CHECK = ."probably"
then
  enable_checks=no
fi

CHECK_LEVEL="0"
if test ".$enable_checks" != ".no"
then
  echo "Run-time checking enabled"
  if test "$enable_checks" = "1"
  then
    echo " -> Level 1 (Basic checking)"
    CXXFLAGS="$CXXFLAGS -DCHECK=1"
    CHECK_LEVEL="1"
  elif test "$enable_checks" = "3"
  then
    CXXFLAGS="$CXXFLAGS -DCHECK=3"
    echo " -> Level 3 (Full checking + stack tracing)"
    CHECK_LEVEL="3"
  else
    echo " -> Level 2 (Basic checking + stack tracing)"
    CXXFLAGS="$CXXFLAGS -DCHECK=2"
    CHECK_LEVEL="2"
  fi
else
  echo "Run-time checking disabled"
fi

if test "$enable_signal" != "no"
then
  echo "Segmentation fault handling enabled"
  CXXFLAGS="$CXXFLAGS -DSIGHANDLE"
else
  echo "Segmentation fault handling disabled"
fi

if test "$enable_color" != "no"
then
  echo "Output coloring enabled"
  CXXFLAGS="$CXXFLAGS -DLOGCOLOR"
else
  echo "Output coloring disabled"
fi

if test ."$enable_track" = ."yes"
then
       CXXFLAGS="$CXXFLAGS -DTRACK"
fi

if test ."$enable_sigfpe" = .yes
then
    CXXFLAGS="$CXXFLAGS -DBOUT_FPE"
fi

echo

BOUT_VERSION=$PACKAGE_VERSION
BOUT_VERSION_MAJOR=${PACKAGE_VERSION%%.*}
BOUT_VERSION_MINOR=${PACKAGE_VERSION#*.}
BOUT_VERSION_DOUBLE=${BOUT_VERSION_MAJOR}.$(echo ${BOUT_VERSION_MINOR}|tr -d .)
CXXFLAGS+=" -DBOUT_VERSION_STRING=\"\\\"$BOUT_VERSION\\\"\" -DBOUT_VERSION_DOUBLE=$BOUT_VERSION_DOUBLE "

#############################################################
# Enable Backtrace if possible
#############################################################

if test ."$enable_backtrace" = .yes || test ."$enable_backtrace" = .maybe
then
    works=yes
    if ! which addr2line &>/dev/null
    then
        AC_MSG_WARN("util addr2line not found")
        works=no
    fi
    AC_CHECK_FUNC(popen, , [works=no])
    AC_TRY_COMPILE([#include <execinfo.h>] , , ,
                   [AC_MSG_WARN("execinfo.h not found")
                    works=no])
    AC_TRY_LINK([extern "C"
                 char backtrace();] ,
                 [return backtrace();], ,
                 [AC_MSG_WARN("backtrace function not found")
                  works=no])
    AC_TRY_COMPILE([#include <dlfcn.h>] , , ,
                   [AC_MSG_WARN("dlfcn.h not found")
                    works=no])
    if test .$works = .yes
    then
        CXXFLAGS="$CXXFLAGS -DBACKTRACE"
    else
        if test ."$enable_backtrace" = .yes
        then
            AC_MSG_ERROR([backtrace requested, but cannot be enabled])
        else
            AC_MSG_WARN("backtrace disabled")
        fi
    fi
fi

#############################################################
# Build into shared object (pic)
#############################################################

if test ."$enable_shared" = .yes
then
    # compile as position independent code.
    # -fpic is apparently faster then -fPIC, but -fPIC works always.
    # From a SO comment (https://stackoverflow.com/a/3544211/3384414):
    #  What's more: I did a little experiment here (on x86_64
    #  platform), -fPIC and -fpic appears to have generated the same
    #  code. It seems they generate a different code only on m68k,
    #  PowerPC and SPARC.
    # Therfore use -fPIC for now
    CXXFLAGS+=" -fPIC"
fi

#############################################################
# Git revision number
#############################################################

rev=`git rev-parse HEAD`
if test "$?" = "0"
then
  dnl Revision not needed for build - is checked at compile time.
  dnl CXXFLAGS="$CXXFLAGS -DREVISION=$rev"
  echo "Revision ID: $rev"
fi

GIT_REVISION=$rev

#############################################################
# FFT routines
#############################################################

echo "Searching for FFTW library"

if type fftw-wisdom > /dev/null 2>&1; then
  echo " -> Found fftw-wisdom executable"
  # Path to fftw-wisdom, then remove fft-wisdom and bin/
  fftw_wisdom=`which fftw-wisdom | xargs -0 dirname | xargs -0 dirname`
    with_fftw="$with_fftw $fftw_wisdom"
fi
BOUT_ADDPATH_CHECK_HEADER(fftw3.h, ,AC_MSG_ERROR('FFTW3 is needed by BOUT++'), $with_fftw)
BOUT_ADDPATH_CHECK_LIB(fftw3, fftw_plan_dft_r2c_1d, ,AC_MSG_ERROR("FFTW3 is needed by BOUT++"), $with_fftw)

echo ""

#############################################################
# netCDF support
#############################################################

NCCONF="" # Configuration script

AS_IF([test "x$with_netcdf" != "xno"],
[
  ##########################################
  # Try to find a valid NetCDF configuration
  #
  # at first, try to find ncconf script

  # Search for NetCDF config scripts, prefer ncxx4-config over nc-config
  # Check if the path to the config script has been supplied directly, otherwise
  # check the path provided by --with-netcdf, appending "/bin" if need
  # be, then check system path
  # Set NCCONF to the full path of the found scropt
  AS_CASE([`basename $with_netcdf 2> /dev/null`],
    ["ncxx4-config"], [NCCONF=$with_netcdf],
    ["nc-config"], [NCCONF=$with_netcdf],
    [AC_PATH_PROGS([NCCONF], [ncxx4-config nc-config], [],
       [$with_netcdf$PATH_SEPARATOR$with_netcdf/bin$PATH_SEPARATOR$PATH])])

  ##########################################
  # Get configuration
  AS_IF([test "x$NCCONF" != "x" ],
  [
     # If we found nc-config rather than ncxx4-config, we need to check if it supports C++
     AS_IF([test `basename $NCCONF` = 'nc-config'],
           [AC_MSG_CHECKING([if $NCCONF has C++4 support])
            nc_has_cpp4=`$NCCONF --has-c++4`
            AC_MSG_RESULT([$nc_has_cpp4])
            AC_MSG_CHECKING([if $NCCONF has C++ support])
            nc_has_cpp=`$NCCONF --has-c++`
            AC_MSG_RESULT([$nc_has_cpp])
           ], [
            nc_has_cpp4="yes"
           ])

     NCINC=`$NCCONF --cflags`
     EXTRA_INCS="$EXTRA_INCS $NCINC"
     AS_IF([test "x$nc_has_cpp4" = "xyes"],
       [
        NCLIB=`$NCCONF --libs`
        CXXFLAGS="$CXXFLAGS -DNCDF4"
        echo " -> NetCDF-4 support enabled"
       ], [
        # nc-config might not *say* it has C++ support, but we can try anyway
        AC_MSG_CHECKING([if we can compile NetCDF with C++])
        # Note netcdf_c++ needed
        NCLIB=`$NCCONF --libs | sed s/-lnetcdf/-lnetcdf_c++\ -lnetcdf/`

        save_LIBS=$LIBS
        save_LDFLAGS=$LDFLAGS
        save_CXXFLAGS=$CXXFLAGS
        AC_LANG_PUSH([C++])
        LIBS="$save_LIBS $NCLIB"
        LDFLAGS="$save_LDFLAGS $NCLIB"
        CXXFLAGS="$save_CXXFLAGS $NCINC"
        AC_LINK_IFELSE(
          [AC_LANG_PROGRAM([
             #include <netcdfcpp.h>
             ], [NcFile file("foo.nc");])],
          [AC_MSG_RESULT([yes])
           echo " -> Legacy NetCDF support enabled"],
          [AC_MSG_RESULT([no])
           AC_MSG_FAILURE([Could not compile NetCDF C++ program!])])
        AC_LANG_POP([C++])
        LIBS=$save_LIBS
        LDFLAGS=$save_LDFLAGS
        CXXFLAGS="$save_CXXFLAGS -DNCDF"
       ])
      EXTRA_LIBS="$EXTRA_LIBS $NCLIB"

      file_formats="$file_formats netCDF"
      NCPATH="found"
      NCFOUND=yes
   ], [
     # if nc-config / ncxx4-config is not found, try to find library directly
     BOUT_MSG_DEBUG([calling bout addpath])
     BOUT_ADDPATH_CHECK_HEADER(netcdfcpp.h,
        BOUT_ADDPATH_CHECK_LIB(netcdf, nc_get_att,
            BOUT_ADDPATH_CHECK_LIB(netcdf_c++, nc_close, NCFOUND=yes, NCFOUND=no, [$with_netcdf]),
            NCFOUND=no, [$with_netcdf]),
        NCFOUND=no, [$with_netcdf])

      AS_IF([test "x$NCFOUND" = "xyes"],
        [
         file_formats="$file_formats netCDF"
         echo " -> Legacy NetCDF support enabled"
         NCPATH="found"
         CXXFLAGS="$CXXFLAGS -DNCDF"
        ], [])
  ])

  AS_IF([test $with_netcdf && test "x$NCFOUND" != "xyes" ], AC_MSG_ERROR([NetCDF requested but not found]), [])
  AS_IF([test "x$NCFOUND" != "xyes"], [AC_MSG_NOTICE([ -> NetCDF support disabled])], [])
], [])

echo ""

#############################################################
# Parallel NetCDF support
#############################################################

PNCPATH=""
if test "$with_pnetcdf" != "no"
then
  echo "Searching for Parallel-NetCDF library"

  if test "$with_pnetcdf" != ""
  then
    # Given a path to the library
    AC_CHECK_FILES($with_pnetcdf/include/pnetcdf.h, PNCPATH=$with_pnetcdf,)
    if test "$PNCPATH" = ""
    then
      echo "parallel-netcdf not found in given directory"
    fi
  fi

  # Find the utilities included with pnetcdf
  if test "$PNCPATH" = ""
  then
    path=`which ncmpidump`
    if test "$path" != ""
    then
      path=`dirname $path`
      path=$path/..
      AC_CHECK_FILES($path/include/pnetcdf.h , PNCPATH=$path,)
    fi
  fi
fi

if test "$PNCPATH" != ""
then
  AC_CHECK_FILES($path/lib/libpnetcdf.a, PNCPATHLIB='lib',
    AC_CHECK_FILES($path/lib64/libpnetcdf.a, PNCPATHLIB='lib64', PNCPATH=''))
fi

if test "$PNCPATH" = ""
then
  echo "Parallel-NetCDF support disabled"
else
  # Set a compile-time flag
  CXXFLAGS="$CXXFLAGS -DPNCDF"
  EXTRA_INCS="$EXTRA_INCS -I$PNCPATH/include"
  EXTRA_LIBS="$EXTRA_LIBS -L$PNCPATH/$PNCPATHLIB -lpnetcdf"

  file_formats="$file_formats Parallel-NetCDF"
  echo "Parallel-NetCDF support enabled"
fi

echo ""

#############################################################
# HDF5 library
#############################################################

HAS_HDF5="no"
AX_LIB_HDF5([serial])
AS_IF([test "$with_hdf5" = "yes"], [
  AC_MSG_NOTICE([Found HDF5])
  CXXFLAGS="$CXXFLAGS -DHDF5"
  EXTRA_INCS="$EXTRA_INCS $HDF5_CPPFLAGS"
  EXTRA_LIBS="$EXTRA_LIBS $HDF5_LDFLAGS $HDF5_LIBS"
  HAS_HDF5="yes"
  file_formats="$file_formats HDF5"
  ], [
  AC_MSG_NOTICE([Serial HDF5 support disabled])
])

HAS_PHDF5="no"
AX_LIB_PARALLELHDF5()
AS_IF([test "$with_parallelhdf5" = "yes"], [
  AC_MSG_NOTICE([Found parallel HDF5])
  CXXFLAGS="$CXXFLAGS -DPHDF5"
  EXTRA_INCS="$EXTRA_INCS $PARALLELHDF5_CPPFLAGS"
  EXTRA_LIBS="$EXTRA_LIBS $PARALLELHDF5_LDFLAGS $PARALLELHDF5_LIBS"
  HAS_PHDF5="yes"
  file_formats="$file_formats Parallel HDF5"
  ], [
  AC_MSG_NOTICE([Parallel HDF5 support disabled])
])

echo ""

#############################################################
# Check file formats
#############################################################

if test "$file_formats" = ""
then
  AC_MSG_ERROR([*** At least one file format must be supported])
else
  echo "Supported file formats:$file_formats"
fi

echo

#############################################################
# LAPACK routines (Used for tri- and band-diagonal solvers)
#############################################################

HAS_LAPACK="no"
if ( ( test "$with_lapack" != "no" ) )
then
    AS_IF( [ test .$with_lapack = .guess || test .$with_lapack = . ],
                 [ lapack_path="" ],
                 [ AS_IF(test .$lapack_path != .yes, [lapack_path=$with_lapack
                                 with_lapack=yes], lapack_path="" )])
    dnl echo "\$with_lapack=$with_lapack"
    dnl echo "\$lapack_pat=$lapack_path"
    BOUT_ADDPATH_CHECK_LIB(blas, zgemm_ , [HAS_BLAS=yes
        BOUT_ADDPATH_CHECK_LIB(lapack, zgbsv_ , HAS_LAPACK=yes
            echo "Using LAPACK"
            CXXFLAGS="$CXXFLAGS -DLAPACK"
            ,
          [ AS_IF( [test .$with_lapack = ".yes"], AC_MSG_ERROR("LAPACK requested but not found.")) ], $lapack_path )
        ]	, [ AS_IF( [test .$with_lapack = ".yes"], AC_MSG_ERROR("LAPACK requested but not found.")) ] , $lapack_path)

    echo ""
fi

#############################################################
# PETSc library
#############################################################

if test "x$with_petsc" != "x" && test "$with_petsc" != "no"
then
### This still needs to check for PETSC_ARCH
  if test "$with_petsc" != "yes"
  then
    PETSC_DIR="$with_petsc"
    AC_CHECK_FILE($PETSC_DIR/conf, PETSC_CONFDIR=${PETSC_DIR}/conf,
      AC_CHECK_FILE($PETSC_DIR/lib/petsc/conf, PETSC_CONFDIR=${PETSC_DIR}/lib/petsc/conf,
        AC_MSG_ERROR([*** --with-petsc=$with_petsc was specified but could not find PETSc distribution])
        ))

  else
    AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/conf, PETSC_CONFDIR=${PETSC_DIR}/conf,
      AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/lib/petsc/conf, PETSC_CONFDIR=${PETSC_DIR}/lib/petsc/conf,
        AC_MSG_ERROR([*** --with-petsc was specified but could not find PETSc distribution])
        ))
  fi
  
  PETSC_VERSION_MAJOR="0"
  PETSC_VERSION_MAJOR=`grep 'define PETSC_VERSION_MAJOR' $PETSC_DIR/include/petscversion.h | awk '{print \$3;}'`
  PETSC_VERSION_MINOR="0"
  PETSC_VERSION_MINOR=`grep 'define PETSC_VERSION_MINOR' $PETSC_DIR/include/petscversion.h | awk '{print \$3;}'`
  PETSC_VERSION_RELEASE="0"
  PETSC_VERSION_RELEASE=`grep 'define PETSC_VERSION_RELEASE' $PETSC_DIR/include/petscversion.h | awk '{print \$3;}'`

  if test $PETSC_VERSION_MAJOR != "3"
  then
    AC_MSG_ERROR([PETSc must be version 3, found $PETSC_VERSION_MAJOR])
  fi

  AC_MSG_CHECKING([PETSc interface version])
  if test "$PETSC_VERSION_MINOR" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_1"
    AC_MSG_RESULT([3.1])
  elif test "$PETSC_VERSION_MINOR" = "2"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_2"
    AC_MSG_RESULT([3.2])
  elif test "$PETSC_VERSION_MINOR" = "3" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_3"
    AC_MSG_RESULT([3.3])
  elif test "$PETSC_VERSION_MINOR" = "4" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_4"
    AC_MSG_RESULT([3.4])
  elif test "$PETSC_VERSION_MINOR" = "5" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_5"
    AC_MSG_RESULT([3.5])
  elif test "$PETSC_VERSION_MINOR" = "5" && test "$PETSC_VERSION_RELEASE" = "0"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_DEV"
    AC_MSG_RESULT([dev])
  elif test "$PETSC_VERSION_MINOR" = "7" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_7"
    AC_MSG_RESULT([3.7 "$PETSC_DIR $PETSC_ARCH $PETSC_CONFDIR"])
  else
    AC_MSG_RESULT([unknown])
    echo "***********************************"
    echo "WARNING: Unrecognised PETSc version"
    echo "MAJOR VERSION: $PETSC_VERSION_MAJOR"
    echo "MINOR VERSION: $PETSC_VERSION_MINOR"
    echo "RELEASE      : $PETSC_VERSION_RELEASE"
    echo ""
    echo "Trying PETSc-dev interface"
    echo "***********************************"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_DEV"
  fi

cat <<MFILE > petscmake$$
PETSC_ARCH = ${PETSC_ARCH}
PETSC_DIR  = ${PETSC_DIR}
include ${PETSC_CONFDIR}/rules
include ${PETSC_CONFDIR}/variables
echo_sundials:
	-@echo \${SUNDIALS_LIB}
MFILE

  AC_MSG_CHECKING([PETSc has SUNDIALS support])
  PETSC_HAS_SUNDIALS=`make -f petscmake$$ echo_sundials | grep -c sundials`

  AS_IF([test "$PETSC_HAS_SUNDIALS" = "0"],
        [PETSC_HAS_SUNDIALS="no"],
        [PETSC_HAS_SUNDIALS="yes"])
  AC_MSG_RESULT([$PETSC_HAS_SUNDIALS])

  rm petscmake$$

  # Set the line to be included in the make.conf file
  PETSC="include ${PETSC_CONFDIR}/variables"
else
  PETSC=
fi
AC_SUBST(PETSC, $PETSC)

if test "$PETSC" != ""
then
  echo "Using PETSC_DIR = $PETSC_DIR"
  echo "Using PETSC_ARCH = $PETSC_ARCH"
  PRECON_SOURCE="$PRECON_SOURCE jstruc.cxx"
  EXTRA_INCS="$EXTRA_INCS \$(PETSC_CC_INCLUDES)"
  EXTRA_LIBS="$EXTRA_LIBS \$(PETSC_LIB)"

  if test "$PETSC_HAS_SUNDIALS" = "yes"
  then
    CXXFLAGS="$CXXFLAGS -DPETSC_HAS_SUNDIALS "
  fi
  echo
fi

#############################################################
# SLEPc library
#############################################################
HAS_SLEPC="no"
SLEPC_VARS=""

if test "$with_slepc" != "" && test "$with_slepc" != "no"
then
  echo "Searching for SLEPc"

  if test "$SLEPC_DIR" = ""
  then
    AC_MSG_ERROR([*** --with-slepc specified but SLEPC_DIR not set])
  fi

  # Check that the config directory exists
  AC_CHECK_FILE($SLEPC_DIR/$PETSC_ARCH/conf, ,
    AC_MSG_ERROR([*** --with-slepc specified but could not find SLEPc distribution])
    )

  echo " -> using SLEPC_DIR=$SLEPC_DIR"

  # Find version number
  AC_CHECK_FILE(${SLEPC_DIR}/include/slepcversion.h, ,
    AC_MSG_ERROR([*** --with-slepc specified but could not read include/slepcversion.h])
    )

  SLEPC_VERSION_MAJOR="0"
  SLEPC_VERSION_MAJOR=`grep 'define SLEPC_VERSION_MAJOR' $SLEPC_DIR/include/slepcversion.h | awk '{print \$3;}'`
  SLEPC_VERSION_MINOR="0"
  SLEPC_VERSION_MINOR=`grep 'define SLEPC_VERSION_MINOR' $SLEPC_DIR/include/slepcversion.h | awk '{print \$3;}'`

  if test "$SLEPC_VERSION_MAJOR" = "3" && test "$SLEPC_VERSION_MINOR" = "4"
  then
    # Version 3.4
    echo " -> Using interface for SLEPc 3.4"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_SLEPC -DBOUT_HAS_SLEPC_3_4"
    HAS_SLEPC="yes"

    EXTRA_INCS="$EXTRA_INCS \$(SLEPC_INCLUDE)"
    EXTRA_LIBS="$EXTRA_LIBS \$(SLEPC_LIB)"
    SLEPC_VARS="include $SLEPC_DIR/conf/slepc_variables"
  else
    AC_MSG_WARN([Error - Unknown SLEPc version])
    echo "  -> No SLEPc support"
    HAS_SLEPC="no"
  fi
  echo
fi
AC_SUBST(SLEPC_VARS, $SLEPC_VARS)

#############################################################
# MUMPS library
#############################################################

if test "$with_mumps" != "" && test "$with_mumps" != "no"
then
  MUMPS=""
  echo "Searching for MUMPS"
  if test "$with_mumps" != "yes"
  then
    # User specified a directory
    MUMPS="$with_mumps"

    # Check for libraries

    MUMPS_INC=$MUMPS/include
    MUMPS_LIB=$MUMPS/lib
  fi

  if test "$MUMPS" = ""
  then
    # Look in some standard locations
    AC_CHECK_FILE(/usr/include/dmumps_c.h, MUMPS="/usr/".
    AC_CHECK_FILE($HOME/local/include/dmumps_c.h, MUMPS="$HOME/local/",))

    # Check for libraries

    MUMPS_INC=$MUMPS/include
    MUMPS_LIB=$MUMPS/lib
  fi

  if test "$MUMPS" = "" && test "$PETSC_DIR" != ""
  then
    # Try PETSc installation

    mumps_dir=`ls $PETSC_DIR/externalpackages/ | grep MUMPS`

    if test "$mumps_dir" != ""
    then
      # PETSc has a MUMPS installation
      MUMPS_INC=$PETSC_DIR/externalpackages/$mumps_dir/include
      MUMPS_LIB=$PETSC_DIR/$PETSC_ARCH/lib
      MUMPS=$PETSC_DIR/externalpackages/$mumps_dir
    fi
  fi

  if test "$MUMPS" != ""
  then
    # Check all files needed
    AC_CHECK_FILES( $MUMPS_INC/dmumps_c.h $MUMPS_LIB/libdmumps.a $MUMPS_LIB/libmumps_common.a $MUMPS_LIB/libpord.a $MUMPS_LIB/libscalapack.a $MUMPS_LIB/libblacs.a $MUMPS_LIB/libparmetis.a $MUMPS_LIB/libmetis.a $MUMPS_LIB/libptscotch.a $MUMPS_LIB/libptesmumps.a, ,[ MUMPS="" ])
  fi

  if test "$MUMPS" = ""
  then
    echo " -> MUMPS not found"
  else
    echo " -> MUMPS support enabled"
    echo "   Includes : $MUMPS_INC"
    echo "   Libraries: $MUMPS_LIB"

    EXTRA_INCS="$EXTRA_INCS -I$MUMPS_INC"
    EXTRA_LIBS="$EXTRA_LIBS -L$MUMPS_LIB -ldmumps -lmumps_common -lpord -lscalapack -lblacs -lparmetis -lmetis -lptscotch -lptesmumps"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_MUMPS"
  fi
  echo
fi

#############################################################
# Solver choice: SUNDIALS' IDA, SUNDIALS' CVODE, PVODE
#############################################################

if ( ( test "$with_sundials" != "" ) && ( test "$with_sundials" != "no" ) )
then
  # Set both IDA and CVODE if not set already

  if test "$with_ida" = ""
  then
    with_ida=$with_sundials
  fi

  if test "$with_cvode" = ""
  then
    with_cvode=$with_sundials
  fi
fi

if ( ( test "$with_ida" != "" ) && ( test "$with_ida" != "no" ) )
then
  echo "Searching for SUNDIALS IDA library"
  if test "$with_ida" = "yes"
  then
    # No path specified. Try using sundials-config
    if type sundials-config > /dev/null 2>&1; then
       IDAINCS=`sundials-config -m ida -t p -l c -s cppflags`
       IDALIBS=`sundials-config -m ida -t p -l c -s libs`
       scr=`which sundials-config`
       echo "=> Using SUNDIALS IDA solver, configured with $scr"
    else
       IDAINCS=""
       IDALIBS=""
    fi
  else
    # Specified with path
    echo "Checking given path to SUNDIALS' IDA solver"
    IDALIBS=$with_ida
    AC_CHECK_FILES( $with_ida/include/ida/ida.h $with_ida/include/ida/ida_spgmr.h $with_ida/include/ida/ida_bbdpre.h $with_ida/include/nvector/nvector_parallel.h $with_ida/include/sundials/sundials_types.h , , [
      IDALIBS=""
    ] )
    AC_CHECK_FILES( $with_ida/lib/libsundials_nvecparallel.a $with_ida/lib/libsundials_ida.a, IDALIBSLIB='lib',
      AC_CHECK_FILES( $with_ida/lib64/libsundials_nvecparallel.a $with_ida/lib64/libsundials_ida.a, IDALIBSLIB='lib64',
        AC_CHECK_FILES( $with_ida/lib/libsundials_nvecparallel.so $with_ida/lib/libsundials_ida.so, IDALIBSLIB='lib',
          AC_CHECK_FILES( $with_ida/lib64/libsundials_nvecparallel.so $with_ida/lib64/libsundials_ida.so, IDALIBSLIB='lib64',
        IDALIBS=""
    ) ) ) )
    if test "$IDALIBS" != ""
    then
    IDAINCS="-I$IDALIBS/include"
    IDALIBS="-L$IDALIBS/$IDALIBSLIB -lsundials_ida -lsundials_nvecparallel"
    echo "=> Using SUNDIALS IDA solver, path $with_ida"
    else
      echo " -> Given path not correct. Finding..."
    fi
  fi

  if test "$IDALIBS" = ""
  then
    # Try some known paths
    AC_CHECK_FILE($HOME/local/include/ida/ida.h, IDALIBS=$HOME/local,
    AC_CHECK_FILE($HOME/include/ida/ida.h, IDALIBS=$HOME,
    AC_CHECK_FILE(/usr/include/ida/ida.h, IDALIBS=/usr,
    AC_CHECK_FILE(/usr/local/include/ida/ida.h, IDALIBS=/usr/local,
    AC_CHECK_FILE(/opt/local/include/ida/ida.h, IDALIBS=/opt/local,
    )))))

    if test "$IDALIBS" != ""
    then
    # Check the rest of the files
    AC_CHECK_FILES( $IDALIBS/include/ida/ida.h $IDALIBS/include/ida/ida_spgmr.h $IDALIBS/include/ida/ida_bbdpre.h $IDALIBS/include/nvector/nvector_parallel.h $IDALIBS/include/sundials/sundials_types.h , , [
          IDALIBS=""
        ] )
    AC_CHECK_FILES( $IDALIBS/lib/libsundials_nvecparallel.a $IDALIBS/lib/libsundials_ida.a, IDALIBSLIB='lib',
      AC_CHECK_FILES( $IDALIBS/lib64/libsundials_nvecparallel.a $IDALIBS/lib64/libsundials_ida.a, IDALIBSLIB='lib64',
        AC_CHECK_FILES( $IDALIBS/lib/libsundials_nvecparallel.so $IDALIBS/lib/libsundials_ida.so, IDALIBSLIB='lib',
          AC_CHECK_FILES( $IDALIBS/lib64/libsundials_nvecparallel.so $IDALIBS/lib64/libsundials_ida.so, IDALIBSLIB='lib64',
            IDALIBS=""
    ) ) ) )
    fi

    if test "$IDALIBS" != ""
    then
    IDAINCS="-I$IDALIBS/include"
    IDALIBS="-L$IDALIBS/$IDALIBSLIB -lsundials_ida -lsundials_nvecparallel"
    fi
  fi

  if test "$IDALIBS" = ""
  then
    echo "=> IDA solver not available"
    echo
  else
    # Compile in the IDA solver
    echo "=> IDA solver enabled"
    echo
    EXTRA_LIBS="$EXTRA_LIBS $IDALIBS"
    EXTRA_INCS="$EXTRA_INCS $IDAINCS"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_IDA"
  fi
fi

if ( ( test "$with_cvode" != "" ) && ( test "$with_cvode" != "no" ) )
then
  echo "Searching for SUNDIALS' CVODE library"
  if test "$with_cvode" = "yes"
  then
    # No path specified. Try using sundials-config
    if type sundials-config > /dev/null 2>&1; then
       CVODEINCS=`sundials-config -m cvode -t p -l c -s cppflags`
       CVODELIBS=`sundials-config -m cvode -t p -l c -s libs`
       scr=`which sundials-config`
       echo "=> Using SUNDIALS CVODE solver, configured with $scr"
    else
       CVODEINCS=""
       CVODELIBS=""
    fi
  else
    # Specified with path
    echo "Checking given path to SUNDIALS' CVODE solver"
    CVODELIBS=$with_cvode
    AC_CHECK_FILES( $with_cvode/include/cvode/cvode.h $with_cvode/include/cvode/cvode_spgmr.h $with_cvode/include/cvode/cvode_bbdpre.h $with_cvode/include/nvector/nvector_parallel.h, , [
      CVODELIBS=""
    ] )
    AC_CHECK_FILES( $with_cvode/lib/libsundials_cvode.a $with_cvode/lib/libsundials_nvecparallel.a, CVODELIBSLIB='lib',
      AC_CHECK_FILES( $with_cvode/lib64/libsundials_cvode.a $with_cvode/lib64/libsundials_nvecparallel.a, CVODELIBSLIB='lib64',
        AC_CHECK_FILES( $with_cvode/lib/libsundials_cvode.so $with_cvode/lib/libsundials_nvecparallel.so, CVODELIBSLIB='lib',
          AC_CHECK_FILES( $with_cvode/lib64/libsundials_cvode.so $with_cvode/lib64/libsundials_nvecparallel.so, CVODELIBSLIB='lib64', [
          CVODELIBS=""
        ]) ) ) )


    if test "$CVODELIBS" != ""
    then
      CVODEINCS="-I$with_cvode/include"
      CVODELIBS="-L$with_cvode/$CVODELIBSLIB -lsundials_cvode -lsundials_nvecparallel"
      echo "=> Using SUNDIALS CVODE solver, path $with_cvode"
    else
      echo " -> Given path not correct. Finding..."
    fi
  fi

  if test "$CVODELIBS" = ""
  then
    # Try some known paths
    AC_CHECK_FILE($HOME/local/include/cvode/cvode_spgmr.h, CVODELIBS=$HOME/local,
    AC_CHECK_FILE($HOME/include/cvode/cvode_spgmr.h, CVODELIBS=$HOME,
    AC_CHECK_FILE(/usr/include/cvode/cvode_spgmr.h, CVODELIBS=/usr,
    AC_CHECK_FILE(/usr/local/include/cvode/cvode_spgmr.h, CVODELIBS=/usr/local,
    AC_CHECK_FILE(/opt/local/include/cvode/cvode_spgmr.h, CVODELIBS=/opt/local,
    )))))

    if test "$CVODELIBS" != ""
    then
    # Check the rest of the files
    AC_CHECK_FILES( $CVODELIBS/include/cvode/cvode.h $CVODELIBS/include/cvode/cvode_bbdpre.h $CVODELIBS/include/nvector/nvector_parallel.h , , [
      CVODELIBS=""
    ])
        AC_CHECK_FILES( $CVODELIBS/lib/libsundials_cvode.a $CVODELIBS/lib/libsundials_nvecparallel.a, CVODELIBSLIB='lib',
          AC_CHECK_FILES( $CVODELIBS/lib64/libsundials_cvode.a $CVODELIBS/lib64/libsundials_nvecparallel.a, CVODELIBSLIB='lib64',
            AC_CHECK_FILES( $CVODELIBS/lib/libsundials_cvode.so $CVODELIBS/lib/libsundials_nvecparallel.so, CVODELIBSLIB='lib',
          AC_CHECK_FILES( $CVODELIBS/lib64/libsundials_cvode.so $CVODELIBS/lib64/libsundials_nvecparallel.so, CVODELIBSLIB='lib64', [
            CVODELIBS=""
          ]) ) ) )
    fi
    if test "$CVODELIBS" != ""
    then
      CVODEINCS="-I$CVODELIBS/include"
      CVODELIBS="-L$CVODELIBS/$CVODELIBSLIB -lsundials_cvode -lsundials_nvecparallel"
    fi
  fi

  if test "$CVODELIBS" = ""
  then
    echo "=> CVODE solver not available"
    echo
  else
    # Compile in the CVODE solver
    echo "=> CVODE solver enabled"
    echo
    EXTRA_LIBS="$EXTRA_LIBS $CVODELIBS"
    EXTRA_INCS="$EXTRA_INCS $CVODEINCS"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_CVODE"


  fi
fi

############## ARKode ####################

if ( ( test "$with_arkode" != "" ) && ( test "$with_arkode" != "no" ) )
then
  echo "Searching for SUNDIALS' ARKODE library"
  if test "$with_arkode" = "yes"
  then
    # No path specified. Try using sundials-config (not currently available in SUNDIALS, but include here for the future)
    if type sundials-config > /dev/null 2>&1; then
       ARKODEINCS=`sundials-config -m arkode -t p -l c -s cppflags`
       ARKODELIBS=`sundials-config -m arkode -t p -l c -s libs`
       scr=`which sundials-config`
       echo "=> Using SUNDIALS ARKODE solver, configured with $scr"
    else
       ARKODEINCS=""
       ARKODELIBS=""
    fi
  else
    # Specified with path
    echo "Checking given path to SUNDIALS' ARKODE solver"
    ARKODELIBS=$with_arkode
    AC_CHECK_FILES( $with_arkode/include/arkode/arkode.h $with_arkode/include/arkode/arkode_spgmr.h $with_arkode/include/arkode/arkode_bbdpre.h $with_arkode/include/nvector/nvector_parallel.h , , [
      ARKODELIBS=""
    ] )
    AC_CHECK_FILES( $ARKODELIBS/lib/libsundials_arkode.a $ARKODELIBS/lib/libsundials_nvecparallel.a, ARKODELIBSLIB='lib',
      AC_CHECK_FILES( $ARKODELIBS/lib64/libsundials_arkode.a $ARKODELIBS/lib64/libsundials_nvecparallel.a, ARKODELIBSLIB='lib64',
        AC_CHECK_FILES( $ARKODELIBS/lib/libsundials_arkode.so $ARKODELIBS/lib/libsundials_nvecparallel.so, ARKODELIBSLIB='lib',
          AC_CHECK_FILES( $ARKODELIBS/lib64/libsundials_arkode.so $ARKODELIBS/lib64/libsundials_nvecparallel.so, ARKODELIBSLIB='lib64',
          [
            ARKODELIBS=""
      ]) ) ) )

    if test "$ARKODELIBS" != ""
    then
      ARKODEINCS="-I$with_arkode/include"
      ARKODELIBS="-L$with_arkode/$ARKODELIBSLIB -lsundials_arkode -lsundials_nvecparallel"
      echo "=> Using SUNDIALS ARKODE solver, path $with_arkode"
    else
      echo " -> Given path not correct. Finding..."
    fi
  fi

  if test "$ARKODELIBS" = ""
  then
    # Try some known paths
    AC_CHECK_FILE($HOME/local/include/arkode/arkode_spgmr.h, ARKODELIBS=$HOME/local,
    AC_CHECK_FILE($HOME/include/arkode/arkode_spgmr.h, ARKODELIBS=$HOME,
    AC_CHECK_FILE(/usr/include/arkode/arkode_spgmr.h, ARKODELIBS=/usr,
    AC_CHECK_FILE(/usr/local/include/arkode/arkode_spgmr.h, ARKODELIBS=/usr/local,
    AC_CHECK_FILE(/opt/local/include/arkode/arkode_spgmr.h, ARKODELIBS=/opt/local,
    )))))

    if test "$ARKODELIBS" != ""
    then
        # Check the rest of the files
        AC_CHECK_FILES( $ARKODELIBS/include/arkode/arkode.h $ARKODELIBS/include/arkode/arkode_bbdpre.h $ARKODELIBS/include/nvector/nvector_parallel.h, , [
          ARKODELIBS=""
        ])
        AC_CHECK_FILES( $ARKODELIBS/lib/libsundials_arkode.a $ARKODELIBS/lib/libsundials_nvecparallel.a, ARKODELIBSLIB='lib',
          AC_CHECK_FILES( $ARKODELIBS/lib64/libsundials_arkode.a $ARKODELIBS/lib64/libsundials_nvecparallel.a, ARKODELIBSLIB='lib64',
            AC_CHECK_FILES( $ARKODELIBS/lib/libsundials_arkode.so $ARKODELIBS/lib/libsundials_nvecparallel.so, ARKODELIBSLIB='lib',
              AC_CHECK_FILES( $ARKODELIBS/lib64/libsundials_arkode.so $ARKODELIBS/lib64/libsundials_nvecparallel.so, ARKODELIBSLIB='lib64',
              [
                ARKODELIBS=""
          ]) ) ) )
    fi

    if test "$ARKODELIBS" != ""
    then
      ARKODEINCS="-I$ARKODELIBS/include"
      ARKODELIBS="-L$ARKODELIBS/$ARKODELIBSLIB -lsundials_arkode -lsundials_nvecparallel"
    fi
  fi
  if test "$ARKODELIBS" = ""
  then
    echo "=> ARKODE solver not available"
    echo
  else
    # Compile in the ARKODE solver
    echo "=> ARKODE solver enabled"
    echo
    EXTRA_LIBS="$EXTRA_LIBS $ARKODELIBS"
    EXTRA_INCS="$EXTRA_INCS $ARKODEINCS"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_ARKODE"


  fi
fi


#############################################################
# Check types for SUNDIALS
#############################################################

AS_IF([test "x$CVODELIBS" != "x"],
[
    BOUT_CHECK_SUNDIALS_TYPE([CVODE],[$CVODEINCS],[CVODEINT],[
#include <cvode/cvode_bbdpre.h>
extern int cvode_bbd_rhs(CVODEINT Nlocal,
               double, N_Vector u, N_Vector du,
               void *user_data);
    ], [CVBBDPrecInit(nullptr, 0, 0, 0, 0, 0, 0, cvode_bbd_rhs, nullptr);])
], [])

AS_IF([test "x$IDALIBS" != "x"],
[
    BOUT_CHECK_SUNDIALS_TYPE([IDA],[$IDAINCS],[IDAINT],[
#include <ida/ida_bbdpre.h>
extern int ida_bbd_rhs(IDAINT Nlocal,
               double, N_Vector u, N_Vector du, N_Vector rr,
               void *user_data);
      ], [IDABBDPrecInit(nullptr, 0, 0, 0, 0, 0, 0, ida_bbd_rhs, nullptr);])
], [])

AS_IF([test "x$ARKODELIBS" != "x"],
[
    BOUT_CHECK_SUNDIALS_TYPE([ARKODE],[$ARKODEINCS],[ARKODEINT],[
#include <arkode/arkode_bbdpre.h>
extern int arkode_bbd_rhs(ARKODEINT Nlocal,
               double, N_Vector u, N_Vector du,
               void *user_data);
      ], [ARKBBDPrecInit(nullptr, 0, 0, 0, 0, 0, 0, arkode_bbd_rhs, nullptr);])
], [])


#############################################################
# Scorep setup
#############################################################

HAS_SCOREP="no"
AS_IF([test "$with_scorep" != "no"], [
  AC_MSG_NOTICE(["Searching for Scorep executable"])

  AS_IF([test "$with_scorep" != "yes"], [
    # Given a path to the library
    AC_CHECK_FILES($with_scorep/scorep, SCOREPPATH=$with_scorep,)
    AS_IF([test "$SCOREPPATH" = ""],[
      AC_MSG_NOTICE(["scorep not found in given directory"])
    ],[])
  ], 
  [])

  # Find the scorep binary
  AS_IF([test "$SCOREPPATH" = ""], [
    path=`which scorep`
    AS_IF([test "$path" != ""], [
      path=`dirname $path`
      path=$path/..
      AC_CHECK_FILES($path/bin/scorep , SCOREPPATH=$path,)
      ], [])
  ],
  [])

  AS_IF([test "$SCOREPPATH" = ""], [
      AC_MSG_NOTICE(["SCOREP support disabled"])
  ],[
      # Set a compile-time flag
      CXXFLAGS="$CXXFLAGS -DBOUT_HAS_SCOREP"
      MPICXX="scorep --user --nocompiler $MPICXX"
      BOUT_CHECK_PRETTYFUNCTION
      HAS_SCOREP="yes"
      AC_MSG_NOTICE(["Scorep support enabled"])
  ])

  ],[
  AC_MSG_NOTICE([Scorep support disabled])
])

echo ""

#############################################################
# Download + Build PVODE '98
#############################################################

AS_MKDIR_P(externalpackages)
AS_MKDIR_P(lib)
AS_MKDIR_P(include)

HAS_PVODE="no"
if ( test "$with_pvode" != "no" )
then
	AS_IF([test "$enable_pvode_openmp" != "no"  ],[
	      AS_IF([test "$enable_openmp" != "no" ],[
    	               PVODE_FLAGS="$CXXFLAGS $OPENMP_CXXFLAGS $CXX11_FLAGS"
		       AC_MSG_NOTICE([PVODE being built with OpenMP support])
		    ],[
		       AC_MSG_ERROR([Cannot enable openmp in PVODE as configuring with OpenMP disabled])
		       ])
	      ],[
	         PVODE_FLAGS="$CXXFLAGS $CXX11_FLAGS"
	         AC_MSG_NOTICE([PVODE being built without OpenMP support])
	])

   # Clean PVODE
   CXX="$MPICXX" CXXFLAGS=$PVODE_FLAGS MKDIR="$MKDIR_P" RANLIB="$RANLIB" $MAKE \
       clean -C externalpackages/PVODE/ >> config-build.log 2>&1

   #make -C externalpackages/PVODE/ clean

    echo "* Building PVODE" | tee -a config-build.log
    echo "*************************************************************" >> config-build.log

    CXX="$MPICXX" CXXFLAGS=$PVODE_FLAGS MKDIR="$MKDIR_P" RANLIB="$RANLIB" \
      AR="ar" ARFLAGS="$ARFLAGS" $MAKE -C externalpackages/PVODE/ >> config-build.log 2>&1

    ex=$?
    if [[ -f externalpackages/PVODE/lib/libpvode.a ]] && test $ex -eq 0
    then
      echo "*************************************************************" >> config-build.log
      echo "* Successfully built PVODE" | tee -a config-build.log
      echo "*	************************************************************" >> config-build.log
      echo "* Installing PVODE into bout++ sourcetree" | tee -a config-build.log

      echo "*************************************************************" | tee -a config-build.log
    else
      AC_MSG_ERROR(Could not build PVODE. See config-build.log for errors)
    fi

  # Set the correct libraries and copy them to bout
  AS_MKDIR_P(include/pvode)
  cp externalpackages/PVODE/precon/pvbbdpre.h externalpackages/PVODE/include/*.h include/pvode
  cp externalpackages/PVODE/lib/*.a lib/
  EXTRA_LIBS="$EXTRA_LIBS -L\$(BOUT_LIB_PATH) -lpvode -lpvpre"
  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PVODE"
  HAS_PVODE="yes"
fi

#############################################################
# Check environment
#############################################################

if ( test "$CXXINCLUDE" != "" )
then
  echo "================================================"
  echo " WARNING: CXXINCLUDE environment variable set to:"
  echo $CXXINCLUDE
  echo " => This will be added to compile commands"
  echo " If this is not intended, then run"
  echo "   export CXXINCLUDE=''"
  echo " before making BOUT++"
  echo "================================================"
fi

#############################################################
# Gather configuration info for bout-config
#############################################################

EXTRA_INCS="${EXTRA_INCS} ${CPPFLAGS}"

PREFIX=$PWD
IDLCONFIGPATH=$PWD/tools/idllib
PYTHONCONFIGPATH=$PWD/tools/pylib

HAS_IDA="yes"
if test "$IDALIBS" = ""
then
  HAS_IDA="no"
fi

HAS_CVODE="yes"
if test "$CVODELIBS" = ""
then
  HAS_CVODE="no"
fi

HAS_ARKODE="yes"
if test "$ARKODELIBS" = ""
then
  HAS_ARKODE="no"
fi

HAS_NETCDF="yes"
if test "$NCPATH" = ""
then
  HAS_NETCDF="no"
fi

HAS_PNETCDF="yes"
if test "$PNCPATH" = ""
then
  HAS_PNETCDF="no"
fi

HAS_MUMPS="yes"
if test "$MUMPS" = ""
then
  HAS_MUMPS="no"
fi

if test "$PETSC" = ""
then
  HAS_PETSC="no"
  PETSC_VERSION=""
  PETSC_RELEASE=""
  PETSC_HAS_SUNDIALS="no"
else
  HAS_PETSC="yes"
  PETSC_VERSION="$PETSC_VERSION_MAJOR.$PETSC_VERSION_MINOR"
  PETSC_RELEASE="$PETSC_VERSION_RELEASE"
fi

#############################################################
# Print configuration info
#############################################################

echo
echo "Configuration summary"  | tee -a config-build.log
echo

if test "$PETSC" = ""
then
  echo "  PETSc support: no" | tee -a config-build.log
else
  echo "  PETSc support: yes (version $PETSC_VERSION_MAJOR.$PETSC_VERSION_MINOR, release = $PETSC_VERSION_RELEASE)" | tee -a config-build.log
  echo "  PETSc has SUNDIALS support: $PETSC_HAS_SUNDIALS" | tee -a config-build.log
fi

echo "  SLEPc support: $HAS_SLEPC" | tee -a config-build.log

echo "  IDA support: $HAS_IDA" | tee -a config-build.log

echo "  CVODE support: $HAS_CVODE" | tee -a config-build.log

echo "  ARKODE support: $HAS_ARKODE" | tee -a config-build.log

echo "  NetCDF support: $HAS_NETCDF" | tee -a config-build.log

echo "  Parallel-NetCDF support: $HAS_PNETCDF" | tee -a config-build.log

echo "  HDF5 support: $HAS_HDF5 (parallel: $HAS_PHDF5)" | tee -a config-build.log

echo "  MUMPS support: $HAS_MUMPS" | tee -a config-build.log

echo "  Lapack support: $HAS_LAPACK" | tee -a config-build.log

echo "  Scorep support: $HAS_SCOREP" | tee -a config-build.log

echo ""
echo
echo "Data analysis configuration"
echo
echo "=== IDL ==="
echo
echo "Make sure that the tools/idllib directory is in your IDL_PATH"
echo "e.g. by adding to your ~/.bashrc file"
echo
echo "    export IDL_PATH=+$PWD/tools/idllib:'<IDL_DEFAULT>':\$IDL_PATH"
echo
echo "=== Python ==="
echo
echo "Make sure that the tools/pylib directory is in your PYTHONPATH"
echo "e.g. by adding to your ~/.bashrc file"
echo
echo "   export PYTHONPATH=$PWD/tools/pylib/:\$PYTHONPATH"
echo
echo "=== Octave ==="
echo
echo "Add the following to ~/.octaverc"
echo "   addpath('$PWD/tools/octave')"
echo
echo "*** Now run '$MAKE' to compile BOUT++ ***"
echo

# Only make.config is altered by configure
AC_OUTPUT(make.config)

#############################################################
# Use a dummy Makefile to get the cflags and ldflags
#
# This is to capture flags from external libraries such
# as PETSc
#############################################################

CONFIG_CFLAGS=`$MAKE cflags -f output.make`
CONFIG_LDFLAGS=`$MAKE ldflags -f output.make`

#############################################################
# Write configuration to bout-config
#############################################################

AC_SUBST(CONFIG_CFLAGS)
AC_SUBST(CONFIG_LDFLAGS)

# Set path to lib and include here.
# If make install is run then that replaces these paths
BOUT_LIB_PATH=$PWD/lib
BOUT_INCLUDE_PATH=$PWD/include
AC_SUBST(BOUT_LIB_PATH)
AC_SUBST(BOUT_INCLUDE_PATH)


AC_SUBST(PREFIX)
AC_SUBST(IDLCONFIGPATH)
AC_SUBST(PYTHONCONFIGPATH)

AC_SUBST(BOUT_VERSION)
AC_SUBST(GIT_REVISION)
AC_SUBST(CHECK_LEVEL)

AC_SUBST(HAS_PVODE)
AC_SUBST(HAS_IDA)
AC_SUBST(HAS_CVODE)
AC_SUBST(HAS_ARKODE)
AC_SUBST(HAS_NETCDF)
AC_SUBST(HAS_PNETCDF)
AC_SUBST(HAS_HDF5)
AC_SUBST(HAS_MUMPS)
AC_SUBST(HAS_LAPACK)
AC_SUBST(HAS_PETSC)
AC_SUBST(HAS_SLEPC)
AC_SUBST(PETSC_VERSION)
AC_SUBST(PETSC_RELEASE)
AC_SUBST(PETSC_HAS_SUNDIALS)
AC_SUBST(HAS_OPENMP)

AC_OUTPUT(bin/bout-config)
chmod a+x bin/bout-config
