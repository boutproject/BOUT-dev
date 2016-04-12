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
# Process this file with autoconf to produce a configure script.
#
# Changelog:
#
# 2010-03-09 Ben Dudson <bd512@york.ac.uk>
#    * Changing to always require FFTW (removing NR routines)
# 2015-08-08 David Schwörer <schword2@mail.dcu.ie>
#    * Searching for libs in lib and lib64
#

BOUT_VERSION="2.0"

AC_PREREQ(2.59)
AC_INIT(BOUT++, BOUT_VERSION, bd512@york.ac.uk)

AC_ARG_WITH(checks, [  --with-checks=no/1/2/3  	   Set run-time checking level])
AC_ARG_WITH(signal, [  --with-signal=no        	   Disable SEGFAULT handling])
AC_ARG_WITH(track,  [  --with-track            	   Enable variable tracking])
AC_ARG_WITH(pdb,    [  --with-pdb              	   Enable support for PDB files])
AC_ARG_WITH(netcdf, [  --with-netcdf           	   Enable support for netCDF files])
AC_ARG_WITH(hdf5,   [  --with-hdf5                 Enable support for HDF5 files])
AC_ARG_WITH(parallelhdf5,[  --with-parallelhdf5    Enable support for parallel HDF5 files])
AC_ARG_WITH(pnetcdf,[  --with-pnetcdf              Set path to Parallel NetCDF library])
AC_ARG_WITH(debug,  [  --with-debug            	   Enable all debugging flags])
AC_ARG_WITH(optimize,[ --with-optimize=default/1/2/3  Enable optimization])
AC_ARG_WITH(ida,    [  --with-ida=/path/to/ida 	   Use the SUNDIALS IDA solver])
AC_ARG_WITH(cvode,  [  --with-cvode            	   Use the SUNDIALS CVODE solver])
AC_ARG_WITH(sundials, [  --with-sundials           Use CVODE and IDA])
AC_ARG_WITH(facets, [  --with-facets           	   Build the interface to FACETS])
AC_ARG_WITH(fftw,   [  --with-fftw             	   Set directory of FFTW3 library])
AC_ARG_WITH(lapack, [  --with-lapack           	   Use the LAPACK library])
AC_ARG_WITH(petsc,  [  --with-petsc            	   Enable PETSc interface])
AC_ARG_WITH(slepc,  [  --with-slepc                Enable SLEPc interface])
AC_ARG_WITH(pvode,  [  --with-pvode	  	   Build and enable PVODE 98 (DEFAULT)])
AC_ARG_WITH(hypre,  [  --with-hypre                Link to the Hypre library])
AC_ARG_WITH(mumps,  [  --with-mumps                Link with MUMPS library for direct matrix inversions])
AC_ARG_WITH(arkode, [  --with-arkode		   Use the SUNDIALS ARKODE solver])

AC_ARG_VAR(EXTRA_INCS,[Extra compile flags])
AC_ARG_VAR(EXTRA_LIBS,[Extra linking flags])

file_formats=""  # Record which file formats are being supported

# Delete the build log from last time
rm config-build.log

AC_ARG_VAR(CXXFLAGS,[Extra compile flags])
AC_ARG_VAR(LDFLAGS,[Extra linking flags])
AC_ARG_VAR(LDLIBS,[Extra linking libraries])

AC_SUBST(MKDIR)
AC_SUBST(EXTRA_INCS)
AC_SUBST(EXTRA_LIBS)

# Adding variables for additional sources
AC_SUBST(PRECON_SOURCE)
AC_SUBST(FACETS_SOURCE)

AC_DEFUN([CREATE_DIR],
        [
          if ! [[ -d $1 ]]
          then
            mkdir $1
          fi
        ])

dnl#############################################################
dnl# Define functions
dnl#############################################################
# Functions inspired by ESPResSo (espressomd.org)
dnl
dnl Try whether the linker / compiler can find the header / libraries
AC_LANG(C++)

AC_DEFUN([BOUT_MSG_DEBUG],[
	dnl Uncomment line below to enable debugging
	dnl AC_MSG_NOTICE([debug: $1])
	])
dnl define the macro to check for libraries
dnl first  argument is the name of the library
dnl second arg is a function of that library
dnl third  arg is to be executed if found
dnl forth  arg is to be executed if not found
dnl fifth  arg is an additional path to check
AC_DEFUN([BOUT_ADDPATH_CHECK_LIB],[
	AC_MSG_CHECKING([for lib$1])
	LDFLAGS_save=$LDFLAGS
	LIBS="$EXTRA_LIBS -l$1"
	BACL_found=no
	AS_IF([test ."$5" = .yes], [extra_prefix=""],[extra_prefix="$5"])
	AC_TRY_LINK([extern "C"
	char $2();] ,[return $2();], [BACL_found=yes ; break;
     		  BOUT_MSG_DEBUG([found $1 without path])
		],)
	if test $BACL_found != yes ; then
		for prefix in $extra_prefix /usr /opt $HOME $HOME/local /usr/local ; do
			for path in $prefix $prefix/lib $prefix/lib64 $prefix/x86_64-linux-gnu ; do
				if test -d $path
				then
					LDFLAGS="-L$path $LDFLAGS_save"
					BOUT_MSG_DEBUG([try link $1 with $path])
					AC_TRY_LINK([extern "C"
						char $2();] ,[return $2();], [BACL_found=yes ; break;
							BOUT_MSG_DEBUG([found $1 with $path])
							],)
				fi
			done
			AS_IF([test .$BACL_found = .yes],break;)
		done
	fi
	if test $BACL_found = yes ; then
		EXTRA_LIBS=$LIBS
		AC_MSG_RESULT(yes)
	else
		LDFLAGS=$LDFLAGS_save
		AC_MSG_RESULT(no)
	fi
	AS_IF([test .$BACL_found = .yes], [$3],[$4])
])


AC_DEFUN([BOUT_ADDPATH_CHECK_HEADER],[
	AC_MSG_CHECKING([for $1])
	CPPFLAGS_save=$CPPFLAGS
	BACH_found=no
	AS_IF([test ."$4" != .yes], [extra_prefix="$4"],[extra_prefix=""])
	AC_TRY_COMPILE([#include <$1>] ,, [BACH_found=yes ; break;
     		BOUT_MSG_DEBUG([found $1 without path])
		],)
	if test $BACH_found != yes ; then
		for prefix in $extra_prefix /usr /opt $HOME $HOME/local /usr/local ; do
			for path in $prefix $prefix/include ; do
				if test -d $path
				then
					CPPFLAGS="$CPPFLAGS_save -I$path"
					BOUT_MSG_DEBUG([try compile $1 with $path])
					AC_TRY_COMPILE([#include <$1>] ,, [BACH_found=yes ; break;
							BOUT_MSG_DEBUG([found $1 with $path])
							],)
				fi
			done
			AS_IF([test .$BACH_found = .yes],break;)
		done
	fi
	if test $BACH_found = yes ; then
		AC_MSG_RESULT(yes)
	else
		CPPFLAGS=$CPPFLAGS_save
		AC_MSG_RESULT(no)
	fi
	AS_IF([test .$BACH_found = .yes], [$2],[$3])
])



#############################################################
# Checks for programs
#############################################################

MKDIR="mkdir -p"
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

#############################################################
# MPI
#############################################################

AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C++], [
  AC_REQUIRE([AC_PROG_CXX])
  AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
  AC_CHECK_PROGS(MPICXX, mpicxx mpiCC mpic++ hcp mpxlC_r mpxlC mpCC cmpic++ CC, $CXX)
  acx_mpi_save_CXX="$CXX"
  CXX="$MPICXX"
  AC_SUBST(MPICXX)
])

if test x = x"$MPILIBS"; then
  AC_LANG_CASE([C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])])
fi

if test x = x"$MPILIBS"; then
  AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
  AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C++], [
  if test x != x"$MPILIBS";
  then
    AC_MSG_CHECKING([for mpi.h])
    AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
      AC_MSG_RESULT(no)])
  fi
])

# AC_LANG_CASE[C++], [CXX="$acx_mpi_save_CXX"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
  $2
  :
else
  ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
  :
fi
])dnl ACX_MPI

AC_LANG([C++])
ACX_MPI()

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_CONST
#AC_TYPE_SIZE_T
#AC_C_VOLATILE

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_VPRINTF
#AC_CHECK_FUNCS([pow sqrt strcasecmp])

# Check for OpenMP support
: ${enable_openmp=no}  # Disable by default
#AC_OPENMP
# Add the flag to both compile and link stages
#CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
#EXTRA_LIBS="$EXTRA_LIBS $OPENMP_CXXFLAGS"

#############################################################
# General Options
#############################################################

OPT_FLAGS=""
if test "$with_debug" != ""
then
  echo "Enabling all debug options"
  with_checks="3"
  # use -Og with available, otherwise fall back to -O0
  OPT_FLAGS="-g -O0 -Og"
elif test $with_optimize ; then
  optimize=$with_optimize
  test $optimize -eq 0 2>/dev/null
  if test $? -eq 2 ; then
    if test $optimize == "default"
    then
      OPT_FLAGS="-O"
    elif test $optimize == "fast"
    then
      OPT_FLAGS="-Ofast -march=native -funroll-loops"
      DISABLE_CHECK=probably
    else
      AC_MSG_ERROR([unrecognized option: --with=optimize=$optimize])
    fi
  else
    if test $optimize -gt 3
    then
      OPT_FLAGS="-Ofast -march=native -funroll-loops"
      DISABLE_CHECK=probably
    elif test $optimize -gt 2
    then
      OPT_FLAGS="-O3 -march=native -funroll-loops"
      DISABLE_CHECK=probably
    elif test $optimize -gt 1
    then
      OPT_FLAGS="-O2 -march=native"
    elif test $optimze -ge 0
    then
      OPT_FLAGS="-O$optimize"
    else
      AC_MSG_ERROR([unrecognized option: --with=optimize=$optimize])
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
if test -z $with_checks && test .$DISABLE_CHECK = ."probably"
then
  with_checks=no
fi

CHECK_LEVEL="0"
if test "$with_checks" != "no"
then
  echo "Run-time checking enabled"
  if test "$with_checks" = "1"
  then
    echo " -> Level 1 (Basic checking)"
    CXXFLAGS="$CXXFLAGS -DCHECK=1"
    CHECK_LEVEL="1"
  elif test "$with_checks" = "3"
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

if test "$with_signal" != "no"
then
  echo "Segmentation fault handling enabled"
  CXXFLAGS="$CXXFLAGS -DSIGHANDLE"
else
  echo "Segmentation fault handling disabled"
fi

if test "$with_track" != ""
then
       CXXFLAGS="$CXXFLAGS -DTRACK"
fi

echo

#############################################################
# Git revision number
#############################################################

rev=`git rev-parse HEAD`
if test "$?" = "0"
then
  # Attach revision info to flags
  CXXFLAGS="$CXXFLAGS -DREVISION=$rev"
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

if test "$with_netcdf" != "no"
then
  ##########################################
  # Try to find a valid NetCDF configuration
  #
  # at first, try to find ncconf script
	
	
  #echo "Searching for NetCDF configuration"

  OLDNCDF="0" # Indicates whether nc-config or ncxx4-config

    if (test "$with_netcdf" = "yes") || (test "$with_netcdf" = "")
    then
      NCCONF=""
      # No path specified. Try using nc-config
      if type ncxx4-config &> /dev/null ; then
        NCCONF=`which ncxx4-config` # New interface
      elif type nc-config &> /dev/null ; then
        # Older interface. May or may not have C++ support
        echo " -> nc-config script found:" `which nc-config`
        if test `nc-config --has-c++` = "yes"
        then
          echo " -> C++ support included"
          NCCONF=`which nc-config`
          OLDNCDF="1"
        else
          echo " -> NO C++ SUPPORT"
        fi
      else
            echo " -> No configure script found. Specify using --with-netcdf option"
      fi
    else
      # see if configure scripts are in path
      if type $with_netcdf/ncxx4-config > /dev/null 2>&1; then
        NCCONF="$with_netcdf/ncxx4-config"
      elif type $with_netcdf/bin/ncxx4-config > /dev/null 2>&1; then
        NCCONF="$with_netcdf/bin/ncxx4-config"
      elif type $with_netcdf/nc-config > /dev/null 2>&1; then
        NCCONF="$with_netcdf/nc-config"
        OLDNCDF="1"
      elif type $with_netcdf/bin/nc-config > /dev/null 2>&1; then
        NCCONF="$with_netcdf/bin/nc-config"
        OLDNCDF="1"
      fi
    fi

    ##########################################
    # Get configuration
    if test "$NCCONF" != ""; then
      echo " -> NetCDF config script found: $NCCONF"

      NCINC=`$NCCONF --cflags`
      EXTRA_INCS="$EXTRA_INCS $NCINC"

      if test "$OLDNCDF" == "1"; then
        NCLIB=`$NCCONF --libs | sed s/-lnetcdf/-lnetcdf_c++\ -lnetcdf/` # Note netcdf_c++ needed
        LDFLAGS="$LDFLAGS -Wl,-rpath,`nc-config --prefix`/lib"

        CXXFLAGS="$CXXFLAGS -DNCDF"
        echo " -> NetCDF support enabled"
      else
        NCLIB=`$NCCONF --libs`

        CXXFLAGS="$CXXFLAGS -DNCDF4"
        echo " -> NetCDF-4 support enabled"
      fi
      EXTRA_LIBS="$EXTRA_LIBS $NCLIB"

      file_formats="$file_formats netCDF"
      NCPATH="found"
      NCFOUND=yes
    fi

  dnl if nc-config / ncxx4-config is not found, try to find library directly
  if test .$NCFOUND != ".yes" ; then
    BOUT_MSG_DEBUG([calling bout addpath])
    BOUT_ADDPATH_CHECK_HEADER(netcdfcpp.h,
      BOUT_ADDPATH_CHECK_LIB(netcdf,nc_get_att,
	BOUT_ADDPATH_CHECK_LIB(netcdf_c++,nc_close, NCFOUND=yes ,NCFOUND=no, \
				$with_netcdf )
      ,NCFOUND=no, $with_netcdf )
    ,NCFOUND=no, $with_netcdf)

  AS_IF( test $NCFOUND = yes, [file_formats="$file_formats netCDF";
    echo " -> NetCDF support enabled";NCPATH="found";
	CXXFLAGS="$CXXFLAGS -DNCDF"]
	,)
  fi

  AS_IF( [test $with_netcdf && test $NCFOUND != yes ] , AC_MSG_ERROR("netcdf requested but not found"),)
  AS_IF( test $NCFOUND != yes, echo " -> NetCDF support disabeld",)

fi

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
    if test "$PNCPATH" == ""
    then
      echo "parallel-netcdf not found in given directory"
    fi
  fi

  # Find the utilities included with pnetcdf
  if test "$PNCPATH" == ""
  then
    path=`which ncmpidump`
    if test "$path" != ""
    then
      path=`dirname $path`
      path=$path/..
      AC_CHECK_FILES($path/include/pnetcdf.h , PNCPATH=$path,)
    fi
  fi

  # Try locate
  if test "$PNCPATH" == ""
  then
    path=`locate pnetcdf.h | head -n 1`
    if test "$path" != ""
    then
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

if test "$PNCPATH" == ""
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
	BOUT_ADDPATH_CHECK_LIB(blas, zgemm_ , HAS_BLAS=yes
		BOUT_ADDPATH_CHECK_LIB(lapack, zgbsv_ , HAS_LAPACK=yes
			echo "Using LAPACK"
			CXXFLAGS="$CXXFLAGS -DLAPACK"
			,
		  [ AS_IF( [test .$with_lapack = ".yes"], AC_MSG_ERROR("LAPACK requested but not found.")) ], $lapack_path )
			, [ AS_IF( [test .$with_lapack = ".yes"], AC_MSG_ERROR("LAPACK requested but not found.")) ] , $lapack_path)
fi

echo ""

#############################################################
# FACETS
#############################################################

if test "$with_facets" != ""
then
    FACETS="$with_facets"
      AC_CHECK_FILE($FACETS/include/FacetsIfc.h, FIFCPATH=$FACETS,
      AC_CHECK_FILE(/usr/include/FacetsIfc.h, FIFCPATH=/usr,
      AC_CHECK_FILE(/usr/local/include/FacetsIfc.h, FIFCPATH=/usr/local,
      AC_CHECK_FILE(/opt/local/include/FacetsIfc.h, FIFCPATH=/opt/local,
      AC_CHECK_FILE($HOME/local/include/FacetsIfc.h, FIFCPATH=$HOME/local,
      AC_CHECK_FILE($HOME/software/facetsifc/include/FacetsIfc.h, FIFCPATH=$HOME/software/facetsifc,
      ))))))
  if test "$FIFCPATH" != ""
  then
    echo "Enabling FACETS interface"
    EXTRA_INCS="$EXTRA_INCS -I$FIFCPATH/include"
    FACETS_SOURCE="$FACETS_SOURCE bout_facets.cxx"
  else
    echo "--with-facets was specified but could not find the FACETS interface"
    exit
  fi
fi

#############################################################
# PETSc library
#############################################################

if test "$with_petsc" != "" && test "$with_petsc" != "no"
then
### This still needs to check for PETSC_ARCH
  if test "$with_petsc" != "yes"
  then
    PETSC_DIR="$with_petsc"

  else
    AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/conf, PETSC_CONFDIR=${PETSC_DIR}/conf,
      AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/lib/petsc/conf, PETSC_CONFDIR=${PETSC_DIR}/lib/petsc/conf,
        echo "Error - can't find PETSc distribution"
        exit))
  fi

  PETSC_VERSION_MAJOR="0"
  PETSC_VERSION_MAJOR=`grep 'define PETSC_VERSION_MAJOR' $PETSC_DIR/include/petscversion.h | awk '{print \$3;}'`
  PETSC_VERSION_MINOR="0"
  PETSC_VERSION_MINOR=`grep 'define PETSC_VERSION_MINOR' $PETSC_DIR/include/petscversion.h | awk '{print \$3;}'`
  PETSC_VERSION_RELEASE="0"
  PETSC_VERSION_RELEASE=`grep 'define PETSC_VERSION_RELEASE' $PETSC_DIR/include/petscversion.h | awk '{print \$3;}'`

  if test $PETSC_VERSION_MAJOR != "3"
  then
    echo "PETSc must be at least version 3"
    exit
  fi

  if test "$PETSC_VERSION_MINOR" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_1"
    echo "Using PETSc 3.1 interface"
  elif test "$PETSC_VERSION_MINOR" = "2"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_2"
    echo "Using PETSc 3.2 interface"
  elif test "$PETSC_VERSION_MINOR" = "3" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_3"
    echo "Using PETSc 3.3 interface"
  elif test "$PETSC_VERSION_MINOR" = "4" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_4"
    echo "Using PETSc 3.4 interface"
  elif test "$PETSC_VERSION_MINOR" = "5" && test "$PETSC_VERSION_RELEASE" = "1"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_3_5"
    echo "Using PETSc 3.5 interface"
  elif test "$PETSC_VERSION_MINOR" = "5" && test "$PETSC_VERSION_RELEASE" = "0"
  then
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_PETSC -DBOUT_HAS_PETSC_DEV"
    echo "Using PETSc-dev interface"
  else
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
echo_cxx:
	-@echo \${CXX}
echo_sundials:
	-@echo \${SUNDIALS_LIB}
MFILE

  CXX="`make -f petscmake$$ echo_cxx`"
  HAS_SUNDIALS=`make -f petscmake$$ echo_sundials | grep -c sundials`
  HAS_CXX=`echo $CXX | grep -ic mp`

  if test "$HAS_SUNDIALS" == "0"
  then
    echo "Warning - PETSc has no SUNDIALS support"
    PETSC_HAS_SUNDIALS="no"
  else
    PETSC_HAS_SUNDIALS="yes"
  fi

  if test "$HAS_CXX" == "0"
  then
    echo "Warning - PETSc must have C++ support"
    #rm petscmake$$
    #exit
  fi


  echo "PETSc used C++ compiler $CXX"
  #MPICXX=$CXX

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

  if test "$PETSC_HAS_SUNDIALS" == "yes"
  then
    CXXFLAGS="$CXXFLAGS -DPETSC_HAS_SUNDIALS "
  fi
fi

echo

#############################################################
# SLEPc library
#############################################################
HAS_SLEPC="no"
SLEPC_VARS=""

if test "$with_slepc" != "" && test "$with_slepc" != "no"
then
  #PETSC_ARCH = ${PETSC_ARCH}
  #SLEPC_DIR = ${SLEPC_DIR}

  echo "Searching for SLEPc"

  if test "$SLEPC_DIR" == ""
  then
    echo "Error - SLEPC_DIR not set"
    exit
  fi

  # Check that the config directory exists
  AC_CHECK_FILE($SLEPC_DIR/$PETSC_ARCH/conf, ,
    echo "Error - can't find SLEPc distribution"
    exit)

  echo " -> using SLEPC_DIR=$SLEPC_DIR"

  # Find version number

  AC_CHECK_FILE(${SLEPC_DIR}/include/slepcversion.h, ,
    echo "Error - can't read include/slepcversion.h"
    exit)

  SLEPC_VERSION_MAJOR="0"
  SLEPC_VERSION_MAJOR=`grep 'define SLEPC_VERSION_MAJOR' $SLEPC_DIR/include/slepcversion.h | awk '{print \$3;}'`
  SLEPC_VERSION_MINOR="0"
  SLEPC_VERSION_MINOR=`grep 'define SLEPC_VERSION_MINOR' $SLEPC_DIR/include/slepcversion.h | awk '{print \$3;}'`

  if test "$SLEPC_VERSION_MAJOR" = "3" && test "$SLEPC_VERSION_MINOR" = "4"
  then
    # Version 3.4
    echo "  -> Using interface for SLEPc 3.4"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_SLEPC -DBOUT_HAS_SLEPC_3_4"
    HAS_SLEPC="yes"

    EXTRA_INCS="$EXTRA_INCS \$(SLEPC_INCLUDE)"
    EXTRA_LIBS="$EXTRA_LIBS \$(SLEPC_LIB)"
    SLEPC_VARS="include $SLEPC_DIR/conf/slepc_variables"
  else
    echo "Error - Unknown SLEPc version"
    echo "  -> No SLEPc support"
    HAS_SLEPC="no"
  fi

fi
AC_SUBST(SLEPC_VARS, $SLEPC_VARS)

echo

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

  if test "$MUMPS" == ""
  then
    # Look in some standard locations
    AC_CHECK_FILE(/usr/include/dmumps_c.h, MUMPS="/usr/".
    AC_CHECK_FILE($HOME/local/include/dmumps_c.h, MUMPS="$HOME/local/",))

    # Check for libraries

    MUMPS_INC=$MUMPS/include
    MUMPS_LIB=$MUMPS/lib
  fi

  if test "$MUMPS" == "" && test "$PETSC_DIR" != ""
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

  if test "$MUMPS" == ""
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
        IDALIBS=""
    ) )
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
           IDALIBS=""
	) )
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
      AC_CHECK_FILES( $with_cvode/lib64/libsundials_cvode.a $with_cvode/lib64/libsundials_nvecparallel.a, CVODELIBSLIB='lib64',  [
      CVODELIBS=""
    ]) )


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
          AC_CHECK_FILES( $CVODELIBS/lib64/libsundials_cvode.a $CVODELIBS/lib64/libsundials_nvecparallel.a, CVODELIBSLIB='lib64',  [
            CVODELIBS=""
          ]) )
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
      AC_CHECK_FILES( $ARKODELIBS/lib64/libsundials_arkode.a $ARKODELIBS/lib64/libsundials_nvecparallel.a, ARKODELIBSLIB='lib64',  [
        ARKODELIBS=""
      ]) )

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
          AC_CHECK_FILES( $ARKODELIBS/lib64/libsundials_arkode.a $ARKODELIBS/lib64/libsundials_nvecparallel.a, ARKODELIBSLIB='lib64',  [
            ARKODELIBS=""
          ]) )
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
# HDF5 library
#############################################################
HDFPATH=""  # Indicates not found
HAS_PHDF5="no"

if test "$with_parallelhdf5" != "no" && test "$with_parallelhdf5" != ""
then
  with_hdf5="$with_parallelhdf5"
fi

if test "$with_hdf5" != "no"
then
  ##########################################
  # Try to find a valid HDF5 configuration

  if (test "$with_hdf5" = "yes") || (test "$with_hdf5" = "")
  then
    # No path specified
    echo "Searching for HDF5 configuration"

    # Try some common HDF5 utilities
    for i in h5cc h5pcc h5diff h5ls
    do
      printf "  -> Looking for $i..."
      if type $i > /dev/null 2>&1; then
        printf "yes\n"
        # Remove 'h5cc' then bin
        HDFPATH=`which $i | xargs -0 dirname | xargs -0 dirname`
        break
      fi
      printf "no\n"
    done

    # Check for include and library files
    AC_CHECK_FILES($HDFPATH/include/hdf5.h ,, [
        echo " -> Can't find HDF5 files"
        HDFPATH=""
      ] )

  else
    HDFPATH=$with_hdf5  # Configure given path

    # Check for include and library files
    AC_CHECK_FILES($HDFPATH/include/hdf5.h ,, [
        echo " ERROR: Given HDF5 path not correct."
        exit
      ] )
  fi


  # Try some known paths
  if test "$HDFPATH" = ""
  then
    AC_CHECK_FILE(/usr/include/hdf5.h, HDFPATH=/usr,
    AC_CHECK_FILE(/usr/local/include/hdf5.h, HDFPATH=/usr/local,
    AC_CHECK_FILE(/opt/local/include/hdf5.h, HDFPATH=/opt/local,
    AC_CHECK_FILE($HOME/local/include/hdf5.h, HDFPATH=$HOME/local,
    AC_CHECK_FILE($HOME/include/hdf5.h, HDFPATH=$HOME,
    )))))
  fi

  if test $HDFPATH
  then
      AC_CHECK_FILES($HDFPATH/lib/libhdf5.a   , HDFLIBPATH='lib',
      AC_CHECK_FILES($HDFPATH/lib64/libhdf5.a , HDFLIBPATH='lib64', HDFPATH=''))
      test $HDFLIBPATH || echo " -> Can't find HDF5 files"
  fi

  # Check if found valid path
  if test "$HDFPATH" != ""
  then
    echo "HDF5 path: $HDFPATH"
    CXXFLAGS="$CXXFLAGS -DHDF5"
    if test "$with_phdf5" != "no" && test "$with_parallelhdf5" != ""
    then
      CXXFLAGS="$CXXFLAGS -DPHDF5"
      HAS_PHDF5="yes"
    fi

    EXTRA_INCS="$EXTRA_INCS -I$HDFPATH/include"
    EXTRA_LIBS="$EXTRA_LIBS -L$HDFPATH/$HDFLIBPATH -lhdf5 -lz"

    file_formats="$file_formats HDF5"
  else
    printf "HDF5 library not found\n\n"
  fi
else
  printf "HDF5 support disabled\n"
fi

echo ""

#############################################################
# PACT library
#############################################################

AC_ARG_VAR(PACT,[Path to the Portable Application Toolkit library])

if test "$with_pdb" = "no"
then
  echo "=> PDB support disabled"
else
  echo "Searching for PDB library"

  if test "$PACT" != ""
  then
    # Path set, check it's correct
    AC_CHECK_FILES( $PACT/include/pdb.h  ,, [
      echo "PACT library location not correct. Finding..."
      PACT=""
    ] )
  else
    echo "PACT variable not set. Searching for library..."
  fi

  if ( test "$with_pdb" != "" ) && ( test "$with_pdb" != "yes" ) && ( test "$with_pdb" != "no" )
  then
    # could be a path
    PACT=$with_pdb
    AC_CHECK_FILES( $PACT/include/pdb.h ,, [
      echo "with-pdb library location not correct. Finding..."
      PACT=""
    ] )
  fi

  ## Try some known paths
  if test "$PACT" = ""
  then
    AC_CHECK_FILE(/usr/include/pdb.h, PACT=/usr/,
    AC_CHECK_FILE(/usr/local/include/pdb.h, PACT=/usr/local/,
    AC_CHECK_FILE(/opt/local/include/pdb.h, PACT=/opt/local/,
    AC_CHECK_FILE($HOME/local/include/pdb.h, PACT=$HOME/local/,
    AC_CHECK_FILE($HOME/local/include/pdb.h, PACT=$HOME/local/,
    )))))


  ## Try to find using locate
  if test "$PACT" = ""
  then
    PACTFILE=`locate pdb.h | head -n 1`
    if test "$PACTFILE" != ""
    then
      PACTFILE=`dirname $PACTFILE`
      PACTFILE=$PACTFILE/..
      AC_CHECK_FILES( $PACTFILE/include/pdb.h ,
        PACT=$PACTFILE,
        PACT="")
    else
      echo "checking locate... no"
    fi
  fi

  if test "$PACT" != ""
  then
    # Check the rest of the PACT libraries exist
    AC_CHECK_FILES( $PACT/lib/libpdb.a   $PACT/lib/libpml.a   $PACT/lib/libscore.a,   PACTPATHLIB='lib',
      AC_CHECK_FILES( $PACT/lib64/libpdb.a $PACT/lib64/libpml.a $PACT/lib64/libscore.a, PACTPATHLIB='lib64', PACT='' ) )
    fi
  fi

  # Test if PACT is set
  if test "$PACT" != ""
  then
    with_pdb="yes"

    # Set include and library paths
    EXTRA_INCS="$EXTRA_INCS -I$PACT/include"
    EXTRA_LIBS="$EXTRA_LIBS -L$PACT/$PACTPATHLIB"
  fi

  ## Add PDB support if available
  if test "$with_pdb" = "yes"
  then
    EXTRA_LIBS="$EXTRA_LIBS -lpdb -lpml -lscore \$(shell $PACT/bin/pact -link)"

    # Set a compile-time flag. PDBF because PDB already used in score.h
    CXXFLAGS="$CXXFLAGS -DPDBF"

    file_formats="$file_formats PDB"
    echo "=> PDB support enabled"
  else
    echo "=> PDB support not available"
  fi
fi

echo

#############################################################
# Check file formats
#############################################################

if test "$file_formats" = ""
then
  echo "ERROR: At least one file format must be supported"
  exit
else
  echo "Supported file formats:$file_formats"
fi

echo

#############################################################
# Hypre
#############################################################

if test "$with_hypre" != "no"
then
  # Try to find a valid Hypre installation

  HYPREPATH=""
  if test "$with_hypre" != "yes"
  then
    # Path specified. Check that needed files exist
    HYPREPATH=$with_hypre
    echo "Looking for HYPRE in $HYPREPATH"
    AC_CHECK_FILES( $HYPREPATH/include/HYPRE.h ,, [
      echo " -> Given path not correct. Finding..."
      HYPREPATH=""
    ] )
    AC_CHECK_FILES( $HYPREPATH/lib/libHYPRE.a ,   HYPREPATHLIB='lib',
    AC_CHECK_FILES( $HYPREPATH/lib64/libHYPRE.a , HYPREPATHLIB='lib64', [
      echo " -> Given path not correct. Finding..."
      HYPREPATH=""
    ] ))
  fi

  if test "$HYPREPATH" = ""
  then
    # try some known paths
    AC_CHECK_FILE(/usr/include/HYPRE.h, HYPREPATH=/usr/,
    AC_CHECK_FILE(/usr/local/include/HYPRE.h, HYPREPATH=/usr/local/,
    AC_CHECK_FILE(/opt/local/include/HYPRE.h, HYPREPATH=/opt/local/,
    AC_CHECK_FILE($HOME/local/include/HYPRE.h, HYPREPATH=$HOME/local,
    AC_CHECK_FILE($HOME/include/HYPRE.h, HYPREPATH=$HOME,
    )))))

    if test "$HYPREPATH" = ""
    then
      # Check for the other files
      AC_CHECK_FILES( $HYPREPATH/lib/libHYPRE.a , HYPREPATHLIB='lib',
        AC_CHECK_FILES( $HYPREPATH/lib64/libHYPRE.a , HYPREPATHLIB='lib64', [
          HYPREPATH=""
        ] ))
    fi
  fi

  if test "$HYPREPATH" != ""
  then
    echo "Support for Hypre enabled"
    EXTRA_INCS="$EXTRA_INCS -I$HYPREPATH/include"
    EXTRA_LIBS="$EXTRA_LIBS -L$HYPREPATH/$HYPREPATHLIB -lHYPRE"

    CFLAGS="$CFLAGS -DBOUT_HAS_HYPRE"
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_HYPRE"
  else
    echo "Support for Hypre disabled"
  fi
  echo
fi


#############################################################
# Check types for SUNDIALS
#############################################################

if test "$CVODELIBS" != ""
then
  echo "checking CVODE types..." | tee -a config-build.log

  cat > test_cvode_compile.cxx <<MYEOF
    #include <cvode/cvode_bbdpre.h>

    #include <stdlib.h>

    extern int cvode_bbd_rhs(CVODEINT Nlocal,
                   double, N_Vector u, N_Vector du,
                   void *user_data);

    int main() {
      CVBBDPrecInit(NULL, 0, 0, 0,
                  0, 0, 0, cvode_bbd_rhs, NULL);
      exit (0);
    }
MYEOF

  # Try to compile
  echo "Trying to compile test_cvode_compile.cxx" >> config-build.log

  $CXX -c test_cvode_compile.cxx $CVODEINCS -DCVODEINT=int  >> config-build.log 2>&1

  # Check return value
  if test "$?" = "0"
  then
    CVODEINT=int
  else
    $CXX -c test_cvode_compile.cxx $CVODEINCS -DCVODEINT=long  >> config-build.log 2>&1
    if test "$?" = "0"
    then
      CVODEINT=long
    else
      echo "HEELP! Cannot compile CVODE with either long or int"
      exit
    fi
  fi
  echo "  CVODE BBD uses '$CVODEINT' type" | tee -a config-build.log

  CXXFLAGS="$CXXFLAGS -DCVODEINT=$CVODEINT"
fi

if test "$IDALIBS" != ""
then
  echo "checking IDA types..." | tee -a config-build.log

  cat > test_ida_compile.cxx <<MYEOF
    #include <ida/ida_bbdpre.h>

    #include <stdlib.h>
    extern int ida_bbd_res(IDAINT Nlocal, double,
		       N_Vector u, N_Vector du, N_Vector rr, void *user_data);

    int main() {
      IDABBDPrecInit(NULL, 0, 0, 0, 0, 0,
			 0, ida_bbd_res, NULL);
      exit (0);
    }
MYEOF

  # Try to compile
  echo "Trying to compile test_ida_compile.cxx" >> config-build.log

  $CXX -c test_ida_compile.cxx $IDAINCS -DIDAINT=int  >> config-build.log 2>&1

  # Check return value
  if test "$?" = "0"
  then
    IDAINT=int
  else
    $CXX -c test_ida_compile.cxx $IDAINCS -DIDAINT=long  >> config-build.log 2>&1
    if test "$?" = "0"
    then
      IDAINT=long
    else
      echo "HEELP! Cannot compile IDA with either long or int"
      exit
    fi
  fi
  echo "  IDA BBD uses '$IDAINT' type" | tee -a config-build.log

  CXXFLAGS="$CXXFLAGS -DIDAINT=$IDAINT"
fi

#### Test ARKode compile ########


if test "$ARKODELIBS" != ""
then
  echo "checking ARKODE types..." | tee -a config-build.log

  cat > test_arkode_compile.cxx <<MYEOF
    #include <arkode/arkode_bbdpre.h>

    #include <stdlib.h>

    extern int arkode_bbd_rhs(ARKODEINT Nlocal,
                   double, N_Vector u, N_Vector du,
                   void *user_data);

    int main() {
      ARKBBDPrecInit(NULL, 0, 0, 0,
                  0, 0, 0, arkode_bbd_rhs, NULL);
      exit (0);
    }
MYEOF

  # Try to compile
  echo "Trying to compile test_arkode_compile.cxx" >> config-build.log

  $CXX -c test_arkode_compile.cxx $ARKODEINCS -DARKODEINT=int  >> config-build.log 2>&1

  # Check return value
  if test "$?" = "0"
  then
    ARKODEINT=int
  else
    $CXX -c test_arkode_compile.cxx $ARKODEINCS -DARKODEINT=long  >> config-build.log 2>&1
    if test "$?" = "0"
    then
      ARKODEINT=long
    else
      echo "HEELP! Cannot compile ARKODE with either long or int"
      exit
    fi
  fi
  echo "  ARKODE BBD uses '$ARKODEINT' type" | tee -a config-build.log

  CXXFLAGS="$CXXFLAGS -DARKODEINT=$ARKODEINT"
fi



#############################################################
# Download + Build PVODE '98
#############################################################

CREATE_DIR(externalpackages)
CREATE_DIR(lib)
CREATE_DIR(include)

HAS_PVODE="no"
if ( test "$with_pvode" != "no" )
then



#    # Get a decent file transfer
#    prog=$(basename "$(which curl)")
#    if test "$prog" == ""
#    then
#      prog=$(basename "$(which wget)")
#    fi
#
#    # This first 'tee' command clears the file out so don't have to delete it if it exists
#    echo "*************************************************************" | tee config-build.log
#    if [[ -f $enable_download_pvode ]]
#    then
#      echo "* Using supplied file $enable_download_pvode" | tee -a config-build.log
#      file=$enable_download_pvode
#    else
#      filename=externalpackages/pvode.tar.gz
#      server=ftp.mcs.anl.gov
#      remotename=/pub/petsc/boutpkgs/pvode.tar.gz
#
#      echo "* Downloading PVODE from $server ... " | tee -a config-build.log
#    if test "$prog" == "curl"
#    then
#      curl -o $filename $server$remotename >> config-build.log 2>&1
#    elif test "$prog" == "wget"
#    then
#      wget -O $filename ftp://$server$remotename >> config-build.log 2>&1
#    else
#      ftp -o $filename -a $servername:$remotename >> config-build.log 2>&1
#    fi
#
#    if test $? != "0"
#    then
#      AC_MSG_ERROR(Error with downloading pvode.tar.gz. Check config-build.log for errors)
#    fi
#    file=externalpackages/pvode.tar.gz
#  fi
#  gunzip -c $file > externalpackages/tmp.tar
#  cd externalpackages
#  tar -x -f tmp.tar >> ../config-build.log 2>&1
#  if test $? != "0"
#  then
#    cd -
#    AC_MSG_ERROR([Cannot extract file: $file])
#  fi
#  cd -

   # Clean PVODE
   CXX="$MPICXX" CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS" MKDIR="$MKDIR" RANLIB="$RANLIB" make clean -C externalpackages/PVODE/precon/ >> config-build.log 2>&1
   CXX="$MPICXX" CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS" MKDIR="$MKDIR" RANLIB="$RANLIB" make clean -C externalpackages/PVODE/source/ >> config-build.log 2>&1

   #make -C externalpackages/PVODE/ clean

    echo "* Building PVODE" | tee -a config-build.log
    echo "*************************************************************" >> config-build.log

    CXX="$MPICXX" CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS" MKDIR="$MKDIR" RANLIB="$RANLIB" make -C externalpackages/PVODE/precon/ >> config-build.log 2>&1
    CXX="$MPICXX" CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS" MKDIR="$MKDIR" RANLIB="$RANLIB" make -C externalpackages/PVODE/source/ >> config-build.log 2>&1

    if [[ -f externalpackages/PVODE/lib/libpvode.a ]]
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
  CREATE_DIR(include/pvode)
  cp externalpackages/PVODE/precon/pvbbdpre.h externalpackages/PVODE/include/*.h include/pvode
  cp externalpackages/PVODE/lib/*.a lib/
  EXTRA_LIBS="$EXTRA_LIBS -L\$(BOUT_TOP)/lib -lpvode -lpvpre"
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

HAS_FACETS="yes"
if test "$FACETS_SOURCE" = ""
then
  HAS_FACETS="no"
fi

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

HAS_HDF5="yes"
if test "$HDFPATH" = ""
then
  HAS_HDF5="no"
fi

HAS_PDB="yes"
if test "$PACT" = ""
then
  HAS_PDB="no"
fi

HAS_HYPRE="yes"
if test "$HYPREPATH" = ""
then
  HAS_HYPRE="no"
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

echo "  FACETS support: $HAS_FACETS"   | tee -a config-build.log

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

echo "  PDB support: $HAS_PDB" | tee -a config-build.log

echo "  Hypre support: $HAS_HYPRE" | tee -a config-build.log

echo "  MUMPS support: $HAS_MUMPS" | tee -a config-build.log

echo ""
echo
echo "Data analysis configuration"
echo
echo "=== IDL ==="
echo
echo "Make sure that the tools/idllib directory is in your IDL_PATH"
echo "e.g. by adding to your ~/.bashrc file"
echo
echo "    export IDL_PATH=$PWD/tools/idllib:\$IDL_PATH"
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

# Only make.config is altered by configure
AC_OUTPUT(make.config)

#############################################################
# Use a dummy Makefile to get the cflags and ldflags
#
# This is to capture flags from external libraries such
# as PETSc
#############################################################

CONFIG_CFLAGS=`make cflags -f output.make`
CONFIG_LDFLAGS=`make ldflags -f output.make`

#############################################################
# Write configuration to bout-config
#############################################################

AC_SUBST(CONFIG_CFLAGS)
AC_SUBST(CONFIG_LDFLAGS)

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
AC_SUBST(HAS_PDB)
AC_SUBST(HAS_HYPRE)
AC_SUBST(HAS_MUMPS)
AC_SUBST(HAS_LAPACK)
AC_SUBST(HAS_PETSC)
AC_SUBST(HAS_SLEPC)
AC_SUBST(PETSC_VERSION)
AC_SUBST(PETSC_RELEASE)
AC_SUBST(PETSC_HAS_SUNDIALS)

AC_OUTPUT(bin/bout-config)
chmod a+x bin/bout-config
