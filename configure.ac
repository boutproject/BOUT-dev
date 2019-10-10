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

AC_PREREQ([2.69])
AC_INIT([BOUT++],[4.2.2],[bd512@york.ac.uk])
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
        [Use the LAPACK library])],,[with_lapack=guess])
AC_ARG_WITH(petsc,        [AS_HELP_STRING([--with-petsc],
        [Enable PETSc interface])],,[with_petsc=no])
AC_ARG_WITH(slepc,        [AS_HELP_STRING([--with-slepc],
        [Enable SLEPc interface])],,[with_slepc=no])
AC_ARG_WITH(pvode,        [AS_HELP_STRING([--with-pvode],
        [Build and enable PVODE 98 (DEFAULT)])],,[])
AC_ARG_WITH(mumps,        [AS_HELP_STRING([--with-mumps],
        [Link with MUMPS library for direct matrix inversions])],,[])
AC_ARG_WITH(arkode,       [AS_HELP_STRING([--with-arkode],
        [Use the SUNDIALS ARKODE solver])],,[])
AC_ARG_WITH(scorep,       [AS_HELP_STRING([--with-scorep],
        [Enable support for scorep based instrumentation])],,[with_scorep=no])

dnl --with-hdf5 flags are set in AX_LIB_{PARALLEL}HDF5

AC_ARG_ENABLE(warnings,     [AS_HELP_STRING([--disable-warnings],
        [Disable compiler warnings])],,[])
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
AC_ARG_ENABLE(optimize,     [AS_HELP_STRING([--enable-optimize=no/1/2/3/4],
        [Enable optimization])],,[])
AC_ARG_ENABLE(sigfpe,       [AS_HELP_STRING([--enable-sigfpe],
        [Enable FloatingPointExceptions])],,[])
AC_ARG_ENABLE(backtrace,    [AS_HELP_STRING([--disable-backtrace],
        [Disable function backtrace])],,[enable_backtrace=maybe])
AC_ARG_ENABLE(shared,       [AS_HELP_STRING([--enable-shared],
        [Enable building bout++ into an shared object])],,[enable_shared=no])
AC_ARG_ENABLE(openmp,       [AS_HELP_STRING([--enable-openmp],
        [Enable building with OpenMP support])],,[enable_openmp=no])
AC_ARG_WITH(openmp_schedule,[AS_HELP_STRING([--with-openmp-schedule=static/dynamic/guided/auto],
        [Set OpenMP schedule (default: static)])],,[with_openmp_schedule=static])
AC_ARG_ENABLE(pvode_openmp, [AS_HELP_STRING([--enable-pvode-openmp],
        [Enable building PVODE with OpenMP support])],,[enable_pvode_openmp=no])

AC_ARG_VAR(EXTRA_INCS,[Extra compile flags])
AC_ARG_VAR(EXTRA_LIBS,[Extra linking flags])

file_formats=""  # Record which file formats are being supported

# Delete the build log from last time
rm -f config-build.log

AC_ARG_VAR(CXXFLAGS,[Extra compile flags])
AC_ARG_VAR(LDFLAGS,[Extra linking flags])
AC_ARG_VAR(LIBS,[Extra linking libraries])
LIBS="$LIBS $LDLIBS"

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

# Autoconf inserts "-g -O2" into flags by default
# Set them to be just "-g", but only if the user hasn't already set CXXFLAGS
# We then put "-O2" back in later, assuming optimisations aren't explicitly disabled
: ${CXXFLAGS="-g"}

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

AS_IF([test "x$enable_openmp" = "xyes"], [
  OPENMP_SCHEDULE=$with_openmp_schedule
  OPENMP_CXXFLAGS="$OPENMP_CXXFLAGS -DOPENMP_SCHEDULE=$OPENMP_SCHEDULE"
])

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
        enable_optimize="no"
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

# Always pass -Werror=unknown-warning-option to get Clang to fail on bad
# flags, otherwise they are always appended to the warn_cxxflags variable,
# and Clang warns on them for every compilation unit.
# If this is passed to GCC, it will explode, so the flag must be enabled
# conditionally.
# This check taken from AX_COMPILER_FLAGS_CXXFLAGS
AX_CHECK_COMPILE_FLAG([-Werror=unknown-warning-option],[
  extra_compiler_flags_test="-Werror=unknown-warning-option"
],[
  extra_compiler_flags_test=""
])


AS_IF([test "x$enable_warnings" != "xno"], [
# Some hopefully sensible default compiler warning flags

  AX_APPEND_COMPILE_FLAGS([ dnl
     -Wall dnl
     -Wextra dnl
     -Wnull-dereference dnl
  ], [CXXFLAGS], [$extra_compiler_flags_test])

# Note we explicitly turn off -Wcast-function-type as PETSc *requires*
# we cast a function to the wrong type in MatFDColoringSetFunction

# Also note that gcc ignores unknown flags of the form "-Wno-warning"
# for backwards compatibility. Therefore we need to add the positive
# form as an additional flag which it will choke on (if it doesn't
# exist). See: https://gcc.gnu.org/wiki/FAQ#wnowarning

  AX_APPEND_COMPILE_FLAGS([ dnl
     -Wno-cast-function-type dnl
  ], [CXXFLAGS], [$extra_compiler_flags_test "-Wcast-function-type"])

], [
  AC_MSG_NOTICE([Compiler warnings disabled])
])

OPT_FLAGS=""
AS_IF([test "$enable_debug" != ""], [
  AC_MSG_NOTICE([Enabling all debug options])
  enable_checks="3"
  # use -Og with available, otherwise fall back to -O0
  OPT_FLAGS="-g -O0 -Og -fno-inline -hipa1"
], [
  AS_IF([test "x$enable_optimize" != "xno"], [
    AS_CASE(["$enable_optimize"],
    ["default" | "yes" | ""],
        [AC_MSG_NOTICE([Enabling default optimisations])
         OPT_FLAGS="-O2"],
    ["fast" | "4"],
        [AC_MSG_NOTICE([Enabling level 4 optimisations])
         OPT_FLAGS="-Ofast -fno-finite-math-only -march=native -funroll-loops"
         DISABLE_CHECK=probably],
    ["3"],
        [AC_MSG_NOTICE([Enabling level 3 optimisations])
         OPT_FLAGS="-O3 -march=native -funroll-loops"
         DISABLE_CHECK=probably],
    ["2"],
        [AC_MSG_NOTICE([Enabling level 2 optimisations])
         OPT_FLAGS="-O2 -march=native"],
    ["1" | "0"],
        [AC_MSG_NOTICE([Enabling level $enable_optimize optimisations])
         OPT_FLAGS="-O$enable_optimize"],
    [
      AC_MSG_ERROR([unrecognized option: --enable-optimize=$enable_optimize])
    ])
  ], [OPT_FLAGS=""])
])

# Append optimisation/debug flags if they work with this compiler
AX_APPEND_COMPILE_FLAGS([ dnl
    $OPT_FLAGS dnl
], [CXXFLAGS], [$extra_compiler_flags_test])

# Disable checks if optimization > 2 is used
AS_IF([test -z $enable_checks && test "x$DISABLE_CHECK" = "xprobably"], [
  enable_checks=no
])

CHECK_LEVEL="0"
AS_IF([test "x$enable_checks" != "xno" && test "x$enable_checks" != "x0"], [
  AC_MSG_NOTICE([Run-time checking enabled])
  AS_CASE([$enable_checks],
    [1], [AC_MSG_NOTICE([ -> Level 1 (Basic checking)])
          CXXFLAGS="$CXXFLAGS -DCHECK=1"
          CHECK_LEVEL="1"],
    [3], [AC_MSG_NOTICE([ -> Level 3 (Full checking + stack tracing)])
          CXXFLAGS="$CXXFLAGS -DCHECK=3"
          CHECK_LEVEL="3"],
    [AC_MSG_NOTICE([ -> Level 2 (Basic checking + stack tracing)])
     CXXFLAGS="$CXXFLAGS -DCHECK=2"
     CHECK_LEVEL="2"])
], [
  AC_MSG_NOTICE([Run-time checking disabled])
])

AS_IF([test "x$enable_signal" != "xno"], [
  AC_MSG_NOTICE([Segmentation fault handling enabled])
  CXXFLAGS="$CXXFLAGS -DSIGHANDLE"
], [
  AC_MSG_NOTICE([Segmentation fault handling disabled])
])

AS_IF([test "x$enable_color" != "xno"], [
  AC_MSG_NOTICE([Output coloring enabled])
  CXXFLAGS="$CXXFLAGS -DLOGCOLOR"
], [
  AC_MSG_NOTICE([Output coloring disabled])
])

AS_IF([test "x$enable_track" = "xyes"], [
  AC_MSG_NOTICE([Field name tracking enabled])
  CXXFLAGS="$CXXFLAGS -DTRACK"
], [
  AC_MSG_NOTICE([Field name tracking disabled])
])

AS_IF([test "x$enable_sigfpe" = "xyes"], [
  AC_MSG_NOTICE([Signaling floating point exceptions enabled])
  CXXFLAGS="$CXXFLAGS -DBOUT_FPE"
], [
  AC_MSG_NOTICE([Signaling floating point exceptions disabled])
])

BOUT_VERSION=$PACKAGE_VERSION
BOUT_VERSION_MAJOR=${PACKAGE_VERSION%%.*}
BOUT_VERSION_MINOR=${PACKAGE_VERSION#*.}
BOUT_VERSION_DOUBLE=${BOUT_VERSION_MAJOR}.$(echo ${BOUT_VERSION_MINOR}|tr -d .)
CXXFLAGS="$CXXFLAGS -DBOUT_VERSION_STRING=\"\\\"$BOUT_VERSION\\\"\" -DBOUT_VERSION_DOUBLE=$BOUT_VERSION_DOUBLE "
#############################################################
# Enable Backtrace if possible
#############################################################

AS_IF([test "x$enable_backtrace" = "xyes" || test "x$enable_backtrace" = "xmaybe"], [
  AC_CHECK_PROG([works], [addr2line], [yes], [no])

  AS_IF([test $works = yes], [
    AC_CHECK_FUNCS([popen backtrace], [works=yes], [works=no; break])
  ])

  AS_IF([test $works = yes], [
    AC_CHECK_HEADERS([execinfo.h dlfcn.h], [works=yes], [works=no; break])
  ])

  AS_IF([test $works = yes], [
    AC_SEARCH_LIBS([dladdr], [dl], [works=yes], [works=no; break])
  ])

  AS_IF([test $works = yes], [
    AC_MSG_NOTICE([Native backtrace enabled])
    CXXFLAGS="$CXXFLAGS -DBACKTRACE"
  ], [
    AS_IF([test "x$enable_backtrace" = "xyes"], [
      AC_MSG_ERROR([backtrace requested, but cannot be enabled])
    ], [
      AC_MSG_WARN([Native backtrace disabled])
    ])
  ])
])

#############################################################
# Build into shared object (pic)
#############################################################

AS_IF([test "x$enable_shared" = "xyes"], [
    # compile as position independent code.
    # -fpic is apparently faster then -fPIC, but -fPIC works always.
    # From a SO comment (https://stackoverflow.com/a/3544211/3384414):
    #  What's more: I did a little experiment here (on x86_64
    #  platform), -fPIC and -fpic appears to have generated the same
    #  code. It seems they generate a different code only on m68k,
    #  PowerPC and SPARC.
    # Therfore use -fPIC for now
    CXXFLAGS="$CXXFLAGS -fPIC"
])

#############################################################
# Git revision number
#############################################################

rev=`git rev-parse HEAD`
AS_IF([test $? = 0], [
  AC_MSG_NOTICE([Git revision: $rev])
  GIT_REVISION=$rev
], [
  GIT_REVISION=
])

#############################################################
# FFT routines
#############################################################

AS_IF([test "x$with_fftw" != "xno"], [
AC_PATH_PROG([fftw_path], [fftw-wisdom], [no], [$with_fftw$PATH_SEPARATOR$PATH])

  AS_IF([test "x$fftw_path" != "xno"], [
    fftw_wisdom0=`AS_DIRNAME(["$fftw_path"])`
    fftw_wisdom=`AS_DIRNAME(["$fftw_wisdom0"])`
    with_fftw="$with_fftw $fftw_wisdom"
  ], AC_MSG_NOTICE([FFTW3 requested but fftw-wisdom not found]))

  BOUT_ADDPATH_CHECK_HEADER(fftw3.h, ,AC_MSG_ERROR([FFTW3 requested but header not found]), $with_fftw)
  BOUT_ADDPATH_CHECK_LIB(fftw3, fftw_plan_dft_r2c_1d, ,AC_MSG_ERROR([FFTW3 requested but library not found]), $with_fftw)

  #If we've reached this point then we're happy that we have FFTW so
  #add our define to CXXFLAGS
  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_FFTW"
  HAS_FFTW="yes"
],
[
AC_MSG_NOTICE([Configuring without FFTW3 is not recommended])
HAS_FFTW="no"
])

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
        AC_MSG_NOTICE([ -> NetCDF-4 support enabled])
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
           AC_MSG_NOTICE([ -> Legacy NetCDF support enabled])],
          [AC_MSG_RESULT([no])
           AC_MSG_FAILURE([*** Could not compile NetCDF C++ program!])])
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
         AC_MSG_NOTICE([ -> Legacy NetCDF support enabled])
         NCPATH="found"
         CXXFLAGS="$CXXFLAGS -DNCDF"
        ], [])
  ])

  AS_IF([test $with_netcdf && test "x$NCFOUND" != "xyes" ], AC_MSG_ERROR([NetCDF requested but not found]), [])
  AS_IF([test "x$NCFOUND" != "xyes"], [AC_MSG_NOTICE([ -> NetCDF support disabled])], [])
], [])

#############################################################
# Parallel NetCDF support
#############################################################

PNCPATH=""
AS_IF([test "$with_pnetcdf" != "no" && test "x$with_pnetcdf" != "x" ], [
  AC_MSG_NOTICE([Searching for Parallel-NetCDF library])

  AS_IF([test "x$with_pnetcdf" != "xyes"], [
    # Given a path to the library
    AC_CHECK_FILES([$with_pnetcdf/include/pnetcdf.h], [PNCPATH=$with_pnetcdf],
      AC_MSG_NOTICE([parallel-netcdf not found in given directory])
    )
  ])

  # Find the utilities included with pnetcdf
  AS_IF([test "x$PNCPATH" = "x"], [
    AS_IF([test "x$with_pnetcdf" = "xyes"], [
      AC_PATH_PROG([NCMPIDUMP_PATH], [ncmpidump])
    ], [
      AC_PATH_PROG([NCMPIDUMP_PATH], [ncmpidump], [$with_pnetcdf$PATH_SEPARATOR$PATH])
    ])
    AS_IF([test "$NCMPIDUMP_PATH" != ""], [
      AC_CHECK_FILES([$NCMPIDUMP_PATH/../include/pnetcdf.h], [PNCPATH=$NCMPIDUMP_PATH/../])
    ])
  ])

  AS_IF([test "x$PNCPATH" != "x"], [
    AC_CHECK_FILES($path/lib/libpnetcdf.a, PNCPATHLIB='lib',
      AC_CHECK_FILES($path/lib64/libpnetcdf.a, PNCPATHLIB='lib64', PNCPATH=''))
  ])

  AS_IF([test "x$PNCPATH" = "x" && test "x$with_pnetcdf" != "x"], [
    AC_MSG_FAILURE([*** Parallel-NetCDF requested but not found])
  ])
])

AS_IF([test "x$PNCPATH" = "x"], [
  AC_MSG_NOTICE([Parallel-NetCDF support disabled])
], [
  # Set a compile-time flag
  CXXFLAGS="$CXXFLAGS -DPNCDF"
  EXTRA_INCS="$EXTRA_INCS -I$PNCPATH/include"
  EXTRA_LIBS="$EXTRA_LIBS -L$PNCPATH/$PNCPATHLIB -lpnetcdf"

  file_formats="$file_formats Parallel-NetCDF"
  AC_MSG_NOTICE([Parallel-NetCDF support enabled])
])

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

#############################################################
# Check file formats
#############################################################

AS_IF([test "x$file_formats" = "x"], [
  AC_MSG_ERROR([*** At least one file format must be supported])
], [
  AC_MSG_NOTICE([Supported file formats:$file_formats])
])

#############################################################
# LAPACK routines (Used for tri- and band-diagonal solvers)
#############################################################

HAS_LAPACK="no"
AS_IF([test "x$with_lapack" != "xno"], [
  AS_IF([test "x$with_lapack" = "xguess" || test "x$with_lapack" = "x" ],
    [lapack_path=""],
    [AS_IF([test "x$with_lapack" != xyes],
       [lapack_path=$with_lapack
        with_lapack=yes],
       [lapack_path=""])
    ])
  BOUT_ADDPATH_CHECK_LIB(blas, zgemm_,
    [HAS_BLAS=yes],
    [AS_IF([test "x$with_lapack" = "xyes"],
       AC_MSG_ERROR([LAPACK requested but couldn't find BLAS]))],
    $lapack_path)
  BOUT_ADDPATH_CHECK_LIB(lapack, zgbsv_,
    [HAS_LAPACK=yes
     AC_MSG_NOTICE([Using LAPACK])
     CXXFLAGS="$CXXFLAGS -DLAPACK"
    ],
    [AS_IF([test "x$with_lapack" = "xyes"],
       AC_MSG_ERROR([LAPACK requested but not found]))],
    $lapack_path)
])

#############################################################
# PETSc library
#############################################################

AS_IF([test "x$with_petsc" != "x" && test "$with_petsc" != "no"], [

# Supplying an argument to "--with-petsc" only works if PETSC_ARCH
# *should* be empty. If it should be non-empty, THIS WILL NOT WORK
  AS_IF([test "$with_petsc" != "yes"], [
    PETSC_DIR="$with_petsc"
    PETSC_ARCH=
  ])

  AC_MSG_NOTICE([Using PETSC_DIR=$PETSC_DIR, PETSC_ARCH=$PETSC_ARCH])

# Define a macro for a nice error message that preserves the
# formatting. Use like:
#
# PETSC_ERROR_MESSAGE
#
# with no identation
  m4_define(PETSC_ERROR_MESSAGE,
  [  You may need to specify PETSC_DIR and PETSC_ARCH like so:
      --with-petsc PETSC_DIR=\$PETSC_DIR PETSC_ARCH=\$PETSC_ARCH
  Also see the help online:
      http://bout-dev.readthedocs.io/en/latest/user_docs/advanced_install.html#petsc
  ])

# PETSc changed the location of the conf directory in 3.5, so we
# need to check both locations
  AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/conf, [
    PETSC_CONFDIR=${PETSC_DIR}/conf
  ], [
    AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/lib/petsc/conf, [
      PETSC_CONFDIR=${PETSC_DIR}/lib/petsc/conf
    ], [
      AC_MSG_ERROR([--with-petsc was specified but could not find PETSc distribution.
PETSC_ERROR_MESSAGE])
    ])
  ])

# We've found an installation, need to check we can use it. First we
# need to be able to check the version number, for which we need
# petscverion.h
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$PETSC_DIR/include"
  AC_CHECK_HEADER([petscversion.h],
    [FOUND_PETSC_HEADER=yes],
    [FOUND_PETSC_HEADER=no]
  )

# This is a terrible hack for some Linux distributions (Fedora) that
# install PETSc in a non-supported fashion. This is really the fault
# of PETSc for their weird method of installation
  AS_IF([test $FOUND_PETSC_HEADER = no], [
    AC_CHECK_FILE([${PETSC_CONFDIR}/petscvariables], [], [
      AC_MSG_ERROR([Unable to find either petscversion.h or petscvariables
PETSC_ERROR_MESSAGE])
    ])

# This relies on the assumption that PETSC_ARCH is empty
    PETSC_CC_INCLUDES=$(grep ^PETSC_CC_INCLUDES ${PETSC_CONFDIR}/petscvariables | cut -d= -f 2-)

    AC_MSG_NOTICE([Looking for petscverion.h using $PETSC_CC_INCLUDES from ${PETSC_CONFDIR}/petscvariables])

    # This is the cache variable set by the previous call to
    # AC_CHECK_HEADER. We need to unset it so we can call the macro
    # again, but now with different CPPFLAGS
    AS_UNSET([ac_cv_header_petscversion_h])

    CPPFLAGS="$CPPFLAGS $PETSC_CC_INCLUDES"
    AC_CHECK_HEADER([petscversion.h], [], [
      AC_MSG_ERROR([Couldn't find or include petscversion.h.
PETSC_ERROR_MESSAGE])
    ])
  ], [])

# Now we have the header we want, we can check the version number
  AC_MSG_CHECKING([PETSc is at least 3.4.0])
  AC_EGREP_CPP([yes], [
    #include <petscversion.h>
    #if PETSC_VERSION_GE(3, 4, 0)
      yes
    #endif
  ], [PETSC_VERSION_OK="yes"],
     [PETSC_VERSION_OK="no"])
  AC_MSG_RESULT([$PETSC_VERSION_OK])
  CPPFLAGS="$save_CPPFLAGS"

  AS_IF([test $PETSC_VERSION_OK = no], [
    AC_MSG_ERROR([PETSc version must be at least 3.4.0])
  ])

# Check if PETSc was compiled with SUNDIALS
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$PETSC_DIR/$PETSC_ARCH/include"
  AC_MSG_CHECKING([PETSc has SUNDIALS support])
  AC_EGREP_CPP([yes], [
    #include <petscconf.h>
    #ifdef PETSC_HAVE_SUNDIALS
      yes
    #endif
  ], [PETSC_HAS_SUNDIALS="yes"],
     [PETSC_HAS_SUNDIALS="no"])
  AC_MSG_RESULT([$PETSC_HAS_SUNDIALS])
  CPPFLAGS="$save_CPPFLAGS"

  AS_IF([test "$PETSC_HAS_SUNDIALS" = "yes"], [
    CXXFLAGS="$CXXFLAGS -DPETSC_HAS_SUNDIALS "
  ])

  # Set the line to be included in the make.conf file
  PETSC_MAKE_INCLUDE="include ${PETSC_CONFDIR}/variables"

  HAS_PETSC="yes"

  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC"
  EXTRA_INCS="$EXTRA_INCS \$(PETSC_CC_INCLUDES)"
  EXTRA_LIBS="$EXTRA_LIBS \$(PETSC_LIB)"

], [
  PETSC_MAKE_INCLUDE=
  HAS_PETSC="no"
  PETSC_HAS_SUNDIALS="no"
])

#############################################################
# SLEPc library
#############################################################

AS_IF([test "x$with_slepc" != "x" && test "$with_slepc" != "no"], [

  AS_IF([test $HAS_PETSC = "no"], [
    AC_MSG_ERROR([--with-slepc specified, but no PETSc detected. Please reconfigure and specify --with-petsc
PETSC_ERROR_MESSAGE])
  ])

# Supplying an argument to "--with-slepc" only works if SLEPC_ARCH
# *should* be empty. If it should be non-empty, THIS WILL NOT WORK
  AS_IF([test "$with_slepc" != "yes"], [
    SLEPC_DIR="$with_slepc"
    SLEPC_ARCH=
  ])

  AC_MSG_NOTICE([Using SLEPC_DIR=$SLEPC_DIR, SLEPC_ARCH=$SLEPC_ARCH])

# Define a macro for a nice error message that preserves the
# formatting. Use like:
#
# SLEPC_ERROR_MESSAGE
#
# with no identation
  m4_define(SLEPC_ERROR_MESSAGE,
  [  You may need to specify SLEPC_DIR and SLEPC_ARCH like so:
      --with-slepc SLEPC_DIR=\$SLEPC_DIR SLEPC_ARCH=\$SLEPC_ARCH
  Also see the help online:
      http://bout-dev.readthedocs.io/en/latest/user_docs/advanced_install.html#slepc
  ])

# Slepc changed the location of the conf directory in 3.5, so we
# need to check both locations
  AC_CHECK_FILE($SLEPC_DIR/$SLEPC_ARCH/conf, [
    SLEPC_CONFDIR=${SLEPC_DIR}/conf
  ], [
    AC_CHECK_FILE($SLEPC_DIR/$SLEPC_ARCH/lib/slepc/conf, [
      SLEPC_CONFDIR=${SLEPC_DIR}/lib/slepc/conf
    ], [
      AC_MSG_ERROR([--with-slepc was specified but could not find Slepc distribution.
SLEPC_ERROR_MESSAGE])
    ])
  ])

# We've found an installation, need to check we can use it. First we
# need to be able to check the version number, for which we need
# slepcverion.h
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$SLEPC_DIR/include"
  AC_CHECK_HEADER([slepcversion.h],
    [FOUND_SLEPC_HEADER=yes],
    [FOUND_SLEPC_HEADER=no]
  )

# This is a terrible hack for some Linux distributions (Fedora) that
# install Slepc in a non-supported fashion. This is really the fault
# of Slepc for their weird method of installation
  AS_IF([test $FOUND_SLEPC_HEADER = no], [
    AC_CHECK_FILE([${SLEPC_CONFDIR}/slepc_variables], [], [
      AC_MSG_ERROR([Unable to find either slepcversion.h or slepc_variables
SLEPC_ERROR_MESSAGE])
    ])

# This relies on the assumption that SLEPC_ARCH is empty
    SLEPC_CC_INCLUDES=$(grep ^SLEPC_CC_INCLUDES ${SLEPC_CONFDIR}/slepc_variables | cut -d= -f 2-)

    AC_MSG_NOTICE([Looking for slepcverion.h using $SLEPC_CC_INCLUDES from ${SLEPC_CONFDIR}/slepc_variables])

    # This is the cache variable set by the previous call to
    # AC_CHECK_HEADER. We need to unset it so we can call the macro
    # again, but now with different CPPFLAGS
    AS_UNSET([ac_cv_header_slepcversion_h])

    CPPFLAGS="$CPPFLAGS $SLEPC_CC_INCLUDES"
    AC_CHECK_HEADER([slepcversion.h], [], [
      AC_MSG_ERROR([Couldn't find or include slepcversion.h.
SLEPC_ERROR_MESSAGE])
    ])
  ], [])

# Now we have the header we want, we can check the version number
  AC_MSG_CHECKING([Slepc is at least 3.4.0])
  AC_EGREP_CPP([yes], [
    #include <slepcversion.h>
    #if SLEPC_VERSION_GE(3, 4, 0)
      yes
    #endif
  ], [SLEPC_VERSION_OK="yes"],
     [SLEPC_VERSION_OK="no"])
  AC_MSG_RESULT([$SLEPC_VERSION_OK])
  CPPFLAGS="$save_CPPFLAGS"

  AS_IF([test $SLEPC_VERSION_OK = no], [
    AC_MSG_ERROR([Slepc version must be at least 3.4.0])
  ])

# Check if Slepc was compiled with SUNDIALS
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$SLEPC_DIR/$SLEPC_ARCH/include"

  # Set the line to be included in the make.conf file
  SLEPC_MAKE_INCLUDE="include ${SLEPC_CONFDIR}/slepc_variables"

  HAS_SLEPC="yes"

  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_SLEPC"
  EXTRA_INCS="$EXTRA_INCS \$(SLEPC_INCLUDE)"
  EXTRA_LIBS="$EXTRA_LIBS \$(SLEPC_LIB)"

], [
  SLEPC_MAKE_INCLUDE=
  HAS_SLEPC="no"
])

#############################################################
# MUMPS library
#############################################################

AS_IF([test "$with_mumps" != "" && test "$with_mumps" != "no"], [
  MUMPS=""
  AC_MSG_NOTICE([Searching for MUMPS])
  AS_IF([test "$with_mumps" != "yes"], [
    # User specified a directory
    MUMPS="$with_mumps"
  ])

  AS_IF([test "$MUMPS" = ""], [
    # Look in some standard locations
    AC_CHECK_FILE(/usr/include/dmumps_c.h, MUMPS="/usr/",
      AC_CHECK_FILE($HOME/local/include/dmumps_c.h, MUMPS="$HOME/local/",))
  ])

  AS_IF([test "$MUMPS" != ""], [
    # Check all files needed
    MUMPS_INC=$MUMPS/include
    MUMPS_LIB=$MUMPS/lib
    AC_CHECK_FILES([$MUMPS_INC/dmumps_c.h $MUMPS_LIB/libdmumps.a
      $MUMPS_LIB/libmumps_common.a $MUMPS_LIB/libpord.a
      $MUMPS_LIB/libscalapack.a
      $MUMPS_LIB/libparmetis.a $MUMPS_LIB/libmetis.a
      $MUMPS_LIB/libptscotch.a $MUMPS_LIB/libptesmumps.a], ,
    [MUMPS=""
     MUMPS_INC=""
     MUMPS_LIB=""
    ])
  ])

  save_LIBS=$LIBS
  save_LDFLAGS=$LDFLAGS
  save_CPPFLAGS=$CPPFLAGS

  AS_IF([test "x$MUMPS_INC" != "x"], [
    LIBS="$save_LIBS -ldmumps -lmumps_common -lpord -lscalapack -lparmetis -lmetis -lptscotch -lptesmumps"
    LDFLAGS="$LDFLAGS -L$MUMPS_LIB"
    CPPFLAGS="$CPPFLAGS -I$MUMPS_INC"
    AC_MSG_NOTICE([Trying the following flags for MUMPS:])
    AC_MSG_NOTICE([   Includes : $CPPFLAGS])
    AC_MSG_NOTICE([   Libraries: $LDFLAGS $LIBS])
  ], [
    AC_MSG_NOTICE([Using user supplied flags for MUMPS])
  ])

  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([if we can compile MUMPS program])
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([
       #include <dmumps_c.h>
       ], [
       DMUMPS_STRUC_C mumps_struc;
       dmumps_c(&mumps_struc);])],
    [AC_MSG_RESULT([yes])
     AC_MSG_NOTICE([ -> MUMPS support enabled])
     MUMPS=yes
     CXXFLAGS="$CXXFLAGS -DBOUT_HAS_MUMPS"],
    [AC_MSG_RESULT([no])
     AC_MSG_FAILURE([*** Could not compile MUMPS program! You may need to supply LIBS, LDFLAGS and CPPFLAGS ])])
  AC_LANG_POP([C++])

  LIBS=$save_LIBS
  LDFLAGS=$save_LDFLAGS
  CPPFLAGS=$save_CPPFLAGS

  AS_IF([test "x$MUMPS_INC" != "x"], [
    AC_MSG_NOTICE([   Includes : $MUMPS_INC])
    AC_MSG_NOTICE([   Libraries: $MUMPS_LIB])


    EXTRA_INCS="$EXTRA_INCS -I$MUMPS_INC"
    EXTRA_LIBS="$EXTRA_LIBS -L$MUMPS_LIB -ldmumps -lmumps_common -lpord -lscalapack -lparmetis -lmetis -lptscotch -lptesmumps"
  ])
])

#############################################################
# Solver choice: SUNDIALS' IDA, SUNDIALS' CVODE, PVODE
#############################################################

AS_IF([test "x$with_sundials" != "x" && test "x$with_sundials" != "xno"], [

  # Now follows a few different checks for the version of SUNDIALS.
  # We need the sundials_config.h header, which comes with all the
  # versions we care about

  # If we've been given a path, look in there first
  AS_IF([test "x$with_sundials" != "xyes"], [
    AC_CHECK_FILE([$with_sundials/include/sundials/sundials_config.h],
      [SUNDIALS_INC=$with_sundials/include], [SUNDIALS_INC=""])
  ])

  # If we've got one, add the include dir to the preprocessor flags
  save_CPPFLAGS=$CPPFLAGS
  AS_IF([test "x$SUNDIALS_INC" != "x"], [CPPFLAGS="-I$SUNDIALS_INC"])

  AC_MSG_CHECKING([for SUNDIALS config header])
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([
      #include "sundials/sundials_config.h"
    ], [])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_FAILURE([*** Could not determine SUNDIALS version])])

  AC_MSG_CHECKING([for SUNDIALS minor version])
  AC_EGREP_CPP([^ *\"? *[12]\.[0-5]\.], [
    #include "sundials/sundials_config.h"
    #ifdef SUNDIALS_PACKAGE_VERSION
    SUNDIALS_PACKAGE_VERSION
    #endif
  ], [sundials_minor_ver="too low"], [sundials_minor_ver=ok])
  AC_MSG_RESULT([$sundials_minor_ver])

  CPPFLAGS=$save_CPPFLAGS

  AS_IF([test "$sundials_minor_ver" = "too low"], [
    AC_MSG_FAILURE([*** Unsupported SUNDIALS version: Requires at least 2.6])
  ])

  # Set both IDA and CVODE if not set already
  AS_IF([test "x$with_ida" = "x"], [
    with_ida=$with_sundials
  ])

  AS_IF([test "x$with_cvode" = "x"], [
    with_cvode=$with_sundials
  ])

  AS_IF([test "x$with_arkode" = "x"], [
    with_arkode=$with_sundials
  ])
])

AS_IF([test "x$with_ida" != "x" && test "x$with_ida" != "xno"], [
  BOUT_FIND_SUNDIALS_MODULE([ida], [
    #include <nvector/nvector_parallel.h>
    #include <ida/ida.h>
    extern void foo(N_Vector);
  ], [IDACreate();])
])

AS_IF([test "x$with_cvode" != "x" && test "x$with_cvode" != "xno"], [
  BOUT_FIND_SUNDIALS_MODULE([cvode], [
    #include <nvector/nvector_parallel.h>
    #include <cvode/cvode.h>
    extern void foo(N_Vector);
  ], [
    #if SUNDIALS_VERSION_MAJOR >= 4
    CVodeCreate(0);
    #else
    CVodeCreate(0, 0);
    #endif
  ])
])

AS_IF([test "x$with_arkode" != "x" && test "x$with_arkode" != "xno"], [
  BOUT_FIND_SUNDIALS_MODULE([arkode], [
    #include <nvector/nvector_parallel.h>
    #if SUNDIALS_VERSION_MAJOR >= 4
    #include <arkode/arkode_arkstep.h>
    #else
    #include <arkode/arkode.h>
    #endif
    extern void foo(N_Vector);
  ], [
    #if SUNDIALS_VERSION_MAJOR >= 4
    ARKStepCreate(0, 0, 0, 0);
    #else
    ARKodeCreate();
    #endif
  ])
])

#############################################################
# Scorep setup
#############################################################

HAS_SCOREP="no"
AS_IF([test "$with_scorep" != "no"], [

  AS_IF([test "$with_scorep" != "yes"], [
    AC_MSG_NOTICE([Searching for Scorep executable in $with_scorep])
    AC_PATH_PROG([SCOREPPATH], [scorep], [], [$with_scorep])
  ], [
    AC_MSG_NOTICE([Searching for Scorep executable])
    AC_PATH_PROG([SCOREPPATH], [scorep], [])
  ])

  AS_IF([test "$SCOREPPATH" = ""], [
      AC_MSG_FAILURE([*** Scorep requested, but executable not found.
Please supply the path using --with-scorep=/path/to/scorep])
  ],[
      # Set a compile-time flag
      CXXFLAGS="$CXXFLAGS -DBOUT_HAS_SCOREP"
      MPICXX="$SCOREPPATH --user --nocompiler $MPICXX"
      BOUT_CHECK_PRETTYFUNCTION
      HAS_SCOREP="yes"
      AC_MSG_NOTICE([Scorep support enabled])
  ])

  ],[
  AC_MSG_NOTICE([Scorep support disabled])
])

#############################################################
# Download + Build PVODE '98
#############################################################

AS_MKDIR_P(externalpackages)
AS_MKDIR_P(lib)
AS_MKDIR_P(include)

HAS_PVODE="no"
AS_IF([test "$with_pvode" != "no"], [
  AS_IF([test "$enable_pvode_openmp" != "no"  ], [
    AS_IF([test "$enable_openmp" != "no" ], [
      PVODE_FLAGS="$CXXFLAGS $OPENMP_CXXFLAGS $CXX11_FLAGS"
      AC_MSG_NOTICE([PVODE being built with OpenMP support])
    ], [
      AC_MSG_ERROR([Cannot enable openmp in PVODE as configuring with OpenMP disabled])
    ])
  ], [
    PVODE_FLAGS="$CXXFLAGS $CXX11_FLAGS"
    AC_MSG_NOTICE([PVODE being built without OpenMP support])
  ])
  # Clean PVODE
  CXX="$MPICXX" CXXFLAGS=$PVODE_FLAGS MKDIR="$MKDIR_P" RANLIB="$RANLIB" $MAKE clean -C externalpackages/PVODE/precon/ >> config-build.log 2>&1
  CXX="$MPICXX" CXXFLAGS=$PVODE_FLAGS MKDIR="$MKDIR_P" RANLIB="$RANLIB" $MAKE clean -C externalpackages/PVODE/source/ >> config-build.log 2>&1

  AC_MSG_NOTICE([Building PVODE])
  echo "* Building PVODE" >> config-build.log
  echo "*************************************************************" >> config-build.log

  CXX="$MPICXX" CXXFLAGS=$PVODE_FLAGS MKDIR="$MKDIR_P" RANLIB="$RANLIB" $MAKE -C externalpackages/PVODE/precon/ >> config-build.log 2>&1
  CXX="$MPICXX" CXXFLAGS=$PVODE_FLAGS MKDIR="$MKDIR_P" RANLIB="$RANLIB" $MAKE -C externalpackages/PVODE/source/ >> config-build.log 2>&1

  AS_IF([test -f externalpackages/PVODE/lib/libpvode.a && test -f externalpackages/PVODE/lib/libpvpre.a], [
    AC_MSG_NOTICE([Successfully built PVODE])
    AC_MSG_NOTICE([Installing PVODE into BOUT++ sourcetree])

    echo "*************************************************************" >> config-build.log
    echo "* Successfully built PVODE" >> config-build.log
    echo "*************************************************************" >> config-build.log
    echo "* Installing PVODE into BOUT++ sourcetree" >> config-build.log
    echo "*************************************************************" >> config-build.log
  ], [
    AC_MSG_ERROR(Could not build PVODE. See config-build.log for errors)
  ])

  # Set the correct libraries and copy them to bout
  AS_MKDIR_P(include/pvode)
  cp -r externalpackages/PVODE/include/pvode include
  cp externalpackages/PVODE/lib/*.a lib/
  EXTRA_LIBS="$EXTRA_LIBS -L\$(BOUT_LIB_PATH) -lpvode -lpvpre"
  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PVODE"
  HAS_PVODE="yes"
])

#############################################################
# Localisation (i18n) with gettext
#############################################################

HAS_NLS="no"
# Use macro to test if Natural Language Support (gettext) is available.
# If available sets:
#  - USE_NLS to "yes"
#  - LIBINTL to the linker options
#  - Modifies CPPFLAGS if needed
AM_GNU_GETTEXT([external])
AS_IF([test "$USE_NLS" = "yes"], [
  AC_MSG_NOTICE([Enabling language support with gettext])
  # Turn the .po files into .mo files
  $MAKE -C locale | tee -a config-build.log 2>&1

  # Note: BOUT_LOCALE_PATH is defined in make.config, and may be changed by `make install`.
  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_GETTEXT -DBOUT_LOCALE_PATH=\$(BOUT_LOCALE_PATH)"

  EXTRA_LIBS="$EXTRA_LIBS $LIBINTL"

  # Set variable substituted into bout-config 
  HAS_NLS="yes"
],[
  AC_MSG_NOTICE([Language support with gettext not available])
]) 

#############################################################
# Check environment
#############################################################

AS_IF([test "$CXXINCLUDE" != ""], [
  AC_MSG_NOTICE([================================================])
  AC_MSG_NOTICE([ WARNING: CXXINCLUDE environment variable set to:])
  AC_MSG_NOTICE([$CXXINCLUDE])
  AC_MSG_NOTICE([ => This will be added to compile commands])
  AC_MSG_NOTICE([ If this is not intended, then run])
  AC_MSG_NOTICE([   export CXXINCLUDE=''])
  AC_MSG_NOTICE([ before making BOUT++])
  AC_MSG_NOTICE([================================================])
])

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

# Only make.config is altered by configure
AC_CONFIG_FILES([make.config])
AC_OUTPUT

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
MPARK_VARIANT_INCLUDE_PATH=$PWD/externalpackages/mpark.variant/include
AC_SUBST(BOUT_LIB_PATH)
AC_SUBST(BOUT_INCLUDE_PATH)
AC_SUBST(MPARK_VARIANT_INCLUDE_PATH)

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
AC_SUBST(PETSC_HAS_SUNDIALS)
AC_SUBST(PETSC_MAKE_INCLUDE)
AC_SUBST(PETSC_DIR)
AC_SUBST(PETSC_ARCH)
AC_SUBST(SLEPC_MAKE_INCLUDE)
AC_SUBST(SLEPC_DIR)
AC_SUBST(SLEPC_ARCH)
AC_SUBST(HAS_OPENMP)
AC_SUBST(HAS_NLS)
AC_SUBST(HAS_FFTW)

AC_CONFIG_FILES([bin/bout-config])
AC_OUTPUT
chmod a+x bin/bout-config

#############################################################
# Print configuration info
#############################################################

AC_MSG_NOTICE([-------------------------])
AC_MSG_NOTICE([  Configuration summary  ])
AC_MSG_NOTICE([-------------------------])

AC_MSG_NOTICE([  PETSc support           : $HAS_PETSC (has SUNDIALS: $PETSC_HAS_SUNDIALS)])
AC_MSG_NOTICE([  SLEPc support           : $HAS_SLEPC])
AC_MSG_NOTICE([  IDA support             : $HAS_IDA])
AC_MSG_NOTICE([  CVODE support           : $HAS_CVODE])
AC_MSG_NOTICE([  ARKODE support          : $HAS_ARKODE])
AC_MSG_NOTICE([  FFTW support            : $HAS_FFTW])
AC_MSG_NOTICE([  NetCDF support          : $HAS_NETCDF])
AC_MSG_NOTICE([  Parallel-NetCDF support : $HAS_PNETCDF])
AC_MSG_NOTICE([  HDF5 support            : $HAS_HDF5 (parallel: $HAS_PHDF5)])
AC_MSG_NOTICE([  MUMPS support           : $HAS_MUMPS])
AC_MSG_NOTICE([  Lapack support          : $HAS_LAPACK])
AC_MSG_NOTICE([  Scorep support          : $HAS_SCOREP])
AC_MSG_NOTICE([  OpenMP support          : $HAS_OPENMP (schedule: $OPENMP_SCHEDULE)])
AC_MSG_NOTICE([  Natural language support: $HAS_NLS (path: $localedir)])

AC_MSG_NOTICE([])
AC_MSG_NOTICE([-------------------------------])
AC_MSG_NOTICE([  Data analysis configuration  ])
AC_MSG_NOTICE([-------------------------------])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([=== IDL ===])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Make sure that the tools/idllib directory is in your IDL_PATH])
AC_MSG_NOTICE([e.g. by adding to your ~/.bashrc file])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([    export IDL_PATH=+$PWD/tools/idllib:'<IDL_DEFAULT>':\$IDL_PATH])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([=== Python ===])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Make sure that the tools/pylib directory is in your PYTHONPATH])
AC_MSG_NOTICE([e.g. by adding to your ~/.bashrc file])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([   export PYTHONPATH=$PWD/tools/pylib/:\$PYTHONPATH])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([*** Now run '$MAKE' to compile BOUT++ ***])
