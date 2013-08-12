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
# 

AC_PREREQ(2.59)
AC_INIT(BOUT++, 1.0, bd512@york.ac.uk)

AC_ARG_WITH(checks, [  --with-checks=no/1/2/3  	   Set run-time checking level])
AC_ARG_WITH(signal, [  --with-signal=no        	   Disable SEGFAULT handling])
AC_ARG_WITH(track,  [  --with-track            	   Enable variable tracking])
AC_ARG_WITH(pdb,    [  --with-pdb              	   Enable support for PDB files])
AC_ARG_WITH(netcdf, [  --with-netcdf           	   Enable support for netCDF files])
AC_ARG_WITH(pnetcdf,[  --with-pnetcdf              Set path to Parallel NetCDF library])
AC_ARG_WITH(debug,  [  --with-debug            	   Enable all debugging flags])
AC_ARG_WITH(ida,    [  --with-ida=/path/to/ida 	   Use the SUNDIALS IDA solver])
AC_ARG_WITH(cvode,  [  --with-cvode            	   Use the SUNDIALS CVODE solver])
AC_ARG_WITH(sundials, [  --with-sundals              Use CVODE and IDA])
AC_ARG_WITH(facets, [  --with-facets           	   Build the interface to FACETS])
AC_ARG_WITH(fftw,   [  --with-fftw             	   Set directory of FFTW3 library])
AC_ARG_WITH(lapack, [  --with-lapack           	   Use the LAPACK library])
AC_ARG_WITH(petsc,  [  --with-petsc            	   Enable PETSc interface])
AC_ARG_WITH(pvode,  [  --with-pvode	  	   Build and enable PVODE 98 (DEFAULT)])
AC_ARG_WITH(hypre,  [  --with-hypre                Link to the Hypre library])
AC_ARG_WITH(mumps,  [  --with-mumps                Link with MUMPS library for direct matrix inversions])

AC_ARG_VAR(EXTRA_INCS,[Extra compile flags])
AC_ARG_VAR(EXTRA_LIBS,[Extra linking flags])

file_formats=""  # Record which file formats are being supported

# Delete the build log from last time
rm config-build.log

AC_ARG_VAR(CXXFLAGS,[Extra compile flags])
AC_ARG_VAR(LDFLAGS,[Extra linking flags])

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

#############################################################
# General Options
#############################################################

if test "$with_debug" != ""
then
  echo "Enabling all debug options"
  with_checks="3"
  CXXFLAGS="$CXXFLAGS -g"
else
  CXXFLAGS="$CXXFLAGS -O"
fi

if test "$with_checks" != "no"
then
  echo "Run-time checking enabled"
  if test "$with_checks" = "1"
  then
    echo " -> Level 1 (Basic checking)"
    CXXFLAGS="$CXXFLAGS -DCHECK=1"
  elif test "$with_checks" = "3"
  then
    CXXFLAGS="$CXXFLAGS -DCHECK=3"
    echo " -> Level 3 (Full checking + stack tracing)"
  else
    echo " -> Level 2 (Basic checking + stack tracing)"
    CXXFLAGS="$CXXFLAGS -DCHECK=2"
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
    PETSC="$with_petsc"
  else
    AC_CHECK_FILE($PETSC_DIR/$PETSC_ARCH/conf, PETSC="include \$(PETSC_DIR)/conf/variables",
    echo "Error - can't find PETSc distribution"
    exit)
  fi

cat <<MFILE > petscmake$$
PETSC_ARCH = ${PETSC_ARCH}
PETSC_DIR  = ${PETSC_DIR}
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
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
  elif test "$PETSC_VERSION_MINOR" = "3" && test "$PETSC_VERSION_RELEASE" = "0"
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
    EXTRA_LIBS="$EXTRA_LIBS -L$MUMPS_LIB -ldmumps -lmumps_common -lpord -lscalapack -lblacs -lblas -lparmetis -lmetis -lptscotch -lptesmumps"
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
    AC_CHECK_FILES( $with_ida/include/ida/ida.h $with_ida/include/ida/ida_spgmr.h $with_ida/include/ida/ida_bbdpre.h $with_ida/include/nvector/nvector_parallel.h $with_ida/include/sundials/sundials_types.h $with_ida/lib/libsundials_ida.a $with_ida/lib/libsundials_nvecparallel.a, , [
      IDALIBS=""
    ] )
    if test "$IDALIBS" != ""
    then
	IDAINCS="-I$IDALIBS/include"
	IDALIBS="-L$IDALIBS/lib -lsundials_ida -lsundials_nvecparallel"
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
	AC_CHECK_FILES( $IDALIBS/include/ida/ida.h $IDALIBS/include/ida/ida_spgmr.h $IDALIBS/include/ida/ida_bbdpre.h $IDALIBS/include/nvector/nvector_parallel.h $IDALIBS/include/sundials/sundials_types.h $IDALIBS/lib/libsundials_ida.a $IDALIBS/lib/libsundials_nvecparallel.a, , [
      	  IDALIBS=""
      	] )
    fi
    
    if test "$IDALIBS" != ""
    then
	IDAINCS="-I$IDALIBS/include"
	IDALIBS="-L$IDALIBS/lib -lsundials_ida -lsundials_nvecparallel"
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
    AC_CHECK_FILES( $with_cvode/include/cvode/cvode.h $with_cvode/include/cvode/cvode_spgmr.h $with_cvode/include/cvode/cvode_bbdpre.h $with_cvode/include/nvector/nvector_parallel.h $with_cvode/lib/libsundials_cvode.a $with_cvode/lib/libsundials_nvecparallel.a, , [
      CVODELIBS=""
    ] )
    
    if test "$CVODELIBS" != ""
    then
      CVODEINCS="-I$with_cvode/include"
      CVODELIBS="-L$with_cvode/lib -lsundials_cvode -lsundials_nvecparallel"
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
	AC_CHECK_FILES( $CVODELIBS/include/cvode/cvode.h $CVODELIBS/include/cvode/cvode_bbdpre.h $CVODELIBS/include/nvector/nvector_parallel.h $CVODELIBS/lib/libsundials_cvode.a $CVODELIBS/lib/libsundials_nvecparallel.a, , [
	  CVODELIBS=""
	])
    fi
    if test "$CVODELIBS" != ""
    then
      CVODEINCS="-I$CVODELIBS/include"
      CVODELIBS="-L$CVODELIBS/lib -lsundials_cvode -lsundials_nvecparallel"
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

#############################################################
# FFT routines
#############################################################

echo "Searching for FFTW library"

if test "$with_fftw" = ""
then
  # No path specified
  FFTPATH=""
  echo " -> Path not given. Searching for library..."
else
  FFTPATH=$with_fftw
  AC_CHECK_FILES( $FFTPATH/include/fftw3.h $FFTPATH/lib/libfftw3.a ,, [
    echo " -> Given path not correct. Finding..."
    FFTPATH=""
  ] )
fi

# Test if FFTW_ROOT is set
if test "$FFTPATH" = ""
then
  if test "$FFTW_ROOT" != ""
  then
    AC_CHECK_FILES( $FFTW_ROOT/include/fftw3.h $FFTW_ROOT/lib/libfftw3.a , FFTPATH=$FFTW_ROOT, )
  fi
fi

# Find fftw-wisdom
if type fftw-wisdom > /dev/null 2>&1; then
  echo " -> Found fftw-wisdom executable"
  # Path to fftw-wisdom, then remove fft-wisdom and bin/
  FFTWISDOM=`which fftw-wisdom | xargs -0 dirname | xargs -0 dirname`
  # Check if the files are there
  AC_CHECK_FILES( $FFTWISDOM/include/fftw3.h $FFTWISDOM/lib/libfftw3.a , FFTPATH=$FFTWISDOM, )
fi

# Try some known paths
if test "$FFTPATH" = ""
then
  AC_CHECK_FILE(/usr/include/fftw3.h, FFTPATH=/usr,
  AC_CHECK_FILE(/usr/local/include/fftw3.h, FFTPATH=/usr/local,
  AC_CHECK_FILE(/opt/local/include/fftw3.h, FFTPATH=/opt/local,
  AC_CHECK_FILE($HOME/local/include/fftw3.h, FFTPATH=$HOME/local,
  ))))

  if test "$FFTPATH" != ""
  then
    # Check for the rest of the files

    # NOTE: AC_CHECK_FILES works here, but AC_CHECK_FILE doesn't (?!)
    AC_CHECK_FILES( $FFTPATH/lib/libfftw3.a , , FFTPATH="" )
    fi
fi

if test "$FFTPATH" = ""
then
  echo "ERROR: FFTW not found. Required by BOUT++"
  exit
else
  echo " -> path $FFTPATH"
  EXTRA_INCS="$EXTRA_INCS -I$FFTPATH/include"
  EXTRA_LIBS="$EXTRA_LIBS -L$FFTPATH/lib -lfftw3"
fi
echo ""

#############################################################
# LAPACK routines (Used for tri- and band-diagonal solvers)
#############################################################

if ( ( test "$with_lapack" != "" ) && ( test "$with_lapack" != "no" ) )
then
  echo "Using LAPACK"
  EXTRA_LIBS="$EXTRA_LIBS -llapack -lblas"
  CXXFLAGS="$CXXFLAGS -DLAPACK"	
fi

#############################################################
# netCDF support
#############################################################

NCCONF="" # Configuration script

if test "$with_netcdf" != "no"
then
  ##########################################
  # Try to find a valid NetCDF configuration 

  echo "Searching for NetCDF configuration"

  OLDNCDF="0" # Indicates whether nc-config or ncxx4-config

  if (test "$with_netcdf" = "yes") || (test "$with_netcdf" = "")
  then
    # No path specified. Try using nc-config
    if type ncxx4-config > /dev/null 2>&1; then
      NCCONF=`which ncxx4-config` # New interface
    elif type nc-config > /dev/null 2>&1; then
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
  else
    # No configuration script. Check the NetCDF path manually
    
    NCPATH=$with_netcdf
    
    if test "$NCPATH" != ""
    then
      # Path set, so check if libraries present
      AC_CHECK_FILES( $NCPATH/include/netcdfcpp.h $NCPATH/lib/libnetcdf_c++.a $NCPATH/lib/libnetcdf.a ,, [
        echo " -> Given path not correct. Finding..."
        NCPATH=""
      ] )
    else
      echo " -> Path not given. Searching for library..."
    fi

    # Try some known paths
    if test "$NCPATH" = ""
    then
      AC_CHECK_FILE(/usr/include/netcdfcpp.h, NCPATH=/usr/,
      AC_CHECK_FILE(/usr/local/include/netcdfcpp.h, NCPATH=/usr/local/,
      AC_CHECK_FILE($HOME/local/include/netcdfcpp.h, NCPATH=$HOME/local,
      )))

      if test "$NCPATH" != ""
      then
        # Check for the rest of the files
        AC_CHECK_FILES( $NCPATH/lib/libnetcdf_c++.a $NCPATH/lib/libnetcdf.a, , NCPATH="")
      fi
    fi

    # Try to find using locate
    if test "$NCPATH" = ""
    then
      NCFILE=`locate netcdfcpp.h | head -n 1`
      if test "$NCFILE" != ""
      then
        NCFILE=`dirname $NCFILE`
        NCFILE=$NCFILE/..
        AC_CHECK_FILES( $NCFILE/include/netcdfcpp.h $NCFILE/lib/libnetcdf_c++.a $NCFILE/lib/libnetcdf.a , NCPATH=$NCFILE, with_netcdf="no")
      else
        echo "checking locate... no"
        with_netcdf="no"
      fi
    fi

    # Finally, if found, add to BOUT++

    if test "$with_netcdf" != "no"
    then
      # Set a compile-time flag
      CXXFLAGS="$CXXFLAGS -DNCDF"

      if test "$NCPATH" != ""; then
        echo " -> path $NCPATH"
        EXTRA_INCS="$EXTRA_INCS -I$NCPATH/include"
        EXTRA_LIBS="$EXTRA_LIBS -L$NCPATH/lib"
        LDFLAGS="$LDFLAGS -Wl,-rpath,$NCPATH/lib"
      fi

      EXTRA_LIBS="$EXTRA_LIBS -lnetcdf_c++ -lnetcdf"

      file_formats="$file_formats netCDF"

      echo " -> NetCDF support enabled"
    else
      echo "NetCDF support disabled"
    fi
  fi

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
    AC_CHECK_FILES($with_pnetcdf/include/pnetcdf.h $with_pnetcdf/lib/libpnetcdf.a, PNCPATH=$with_pnetcdf,)
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
      AC_CHECK_FILES($path/include/pnetcdf.h $path/lib/libpnetcdf.a, PNCPATH=$path,)
    fi
  fi

  # Try locate
  if test "$PNCPATH" == ""
  then
    path=`locate pnetcdf.h | head -n 1`
    if test "$path" != ""
    then
      path=$path/..
      AC_CHECK_FILES($path/include/pnetcdf.h $path/lib/libpnetcdf.a, PNCPATH=$path,)
    fi
  fi
fi

if test "$PNCPATH" == ""
then
  echo "Parallel-NetCDF support disabled"
else
  # Set a compile-time flag
  CXXFLAGS="$CXXFLAGS -DPNCDF"
  EXTRA_INCS="$EXTRA_INCS -I$PNCPATH/include"
  EXTRA_LIBS="$EXTRA_LIBS -L$PNCPATH/lib -lpnetcdf"
  
  file_formats="$file_formats Parallel-NetCDF"
  echo "Parallel-NetCDF support enabled"
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
    AC_CHECK_FILES( $PACT/include/pdb.h $PACT/lib/libpdb.a $PACT/lib/libpml.a $PACT/lib/libscore.a ,, [
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

    AC_CHECK_FILES( $PACT/include/pdb.h $PACT/lib/libpdb.a $PACT/lib/libpml.a $PACT/lib/libscore.a ,, [
      echo "with-pdb library location not correct. Finding..."
      PACT=""
    ] )
  fi

  ## Try some known paths
  if test "$PACT" = ""
  then
    AC_CHECK_FILE(/usr/include/pdb.h, PACT=/usr/,
    AC_CHECK_FILE(/usr/local/include/pdb.h, PACT=/usr/local/,
    AC_CHECK_FILE($HOME/local/include/pdb.h, PACT=$HOME/local/,
    AC_CHECK_FILE($HOME/local/include/pdb.h, PACT=$HOME/local/,
    ))))

  if test "$PACT" != ""
  then
    # Check the rest of the PACT libraries exist
    AC_CHECK_FILES( $PACT/lib/libpdb.a $PACT/lib/libpml.a $PACT/lib/libscore.a, , PACT="")
    fi
  fi

  ## Try to find using locate
  if test "$PACT" = ""
  then
    PACTFILE=`locate pdb.h | head -n 1`
    if test "$PACTFILE" != ""
    then
      PACTFILE=`dirname $PACTFILE`
      PACTFILE=$PACTFILE/..
      AC_CHECK_FILES( $PACTFILE/include/pdb.h $PACTFILE/lib/libpdb.a $PACTFILE/lib/libpml.a $PACTFILE/lib/libscore.a,
        PACT=$PACTFILE,
        PACT="")
    else
      echo "checking locate... no"
    fi
  fi

  # Test if PACT is set
  if test "$PACT" != ""
  then
    with_pdb="yes"

    # Set include and library paths
    EXTRA_INCS="$EXTRA_INCS -I$PACT/include"
    EXTRA_LIBS="$EXTRA_LIBS -L$PACT/lib"
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
    AC_CHECK_FILES( $HYPREPATH/include/HYPRE.h $HYPREPATH/lib/libHYPRE.a ,, [
      echo " -> Given path not correct. Finding..."
      HYPREPATH=""
    ] )
  fi

  if test "$HYPREPATH" = ""
  then
    # try some known paths
    AC_CHECK_FILE(/usr/include/HYPRE.h, HYPREPATH=/usr/,
    AC_CHECK_FILE(/usr/local/include/HYPRE.h, HYPREPATH=/usr/local/,
    AC_CHECK_FILE(/opt/local/include/HYPRE.h, HYPREPATH=/opt/local/,
    AC_CHECK_FILE($HOME/local/include/HYPRE.h, HYPREPATH=$HOME/local,
    ))))
    
    if test "$HYPREPATH" = ""
    then
      # Check for the other files
      AC_CHECK_FILES( $HYPREPATH/lib/libHYPRE.a ,, [
        HYPREPATH=""
        ] )
    fi
  fi
  
  if test "$HYPREPATH" != ""
  then
    echo "Support for Hypre enabled"
    EXTRA_INCS="$EXTRA_INCS -I$HYPREPATH/include"
    EXTRA_LIBS="$EXTRA_LIBS -L$HYPREPATH/lib -lHYPRE"
    
    CFLAGS="$CFLAGS -DBOUT_HAS_HYPRE"
  else
    echo "Support for Hypre disabled"
  fi
  echo ""
fi

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
    AC_CHECK_FILES( $HYPREPATH/include/HYPRE.h $HYPREPATH/lib/libHYPRE.a ,, [
      echo " -> Given path not correct. Finding..."
      HYPREPATH=""
    ] )
  fi

  if test "$HYPREPATH" = ""
  then
    # try some known paths
    AC_CHECK_FILE(/usr/include/HYPRE.h, HYPREPATH=/usr/,
    AC_CHECK_FILE(/usr/local/include/HYPRE.h, HYPREPATH=/usr/local/,
    AC_CHECK_FILE(/opt/local/include/HYPRE.h, HYPREPATH=/opt/local/,
    AC_CHECK_FILE($HOME/local/include/HYPRE.h, HYPREPATH=$HOME/local,
    ))))
    
    if test "$HYPREPATH" = ""
    then
      # Check for the other files
      AC_CHECK_FILES( $HYPREPATH/lib/libHYPRE.a ,, [
        HYPREPATH=""
        ] )
    fi
  fi
  
  if test "$HYPREPATH" != ""
  then
    echo "Support for Hypre enabled"
    EXTRA_INCS="$EXTRA_INCS -I$HYPREPATH/include"
    EXTRA_LIBS="$EXTRA_LIBS -L$HYPREPATH/lib -lHYPRE"
    
    CXXFLAGS="$CXXFLAGS -DBOUT_HAS_HYPRE"
  else
    echo "Support for Hypre disabled"
  fi
  echo ""
fi

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
    AC_TRY_COMPILE([#include &lt;mpi.h&gt;],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
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
AC_OPENMP
# Add the flag to both compile and link stages
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
EXTRA_LIBS="$EXTRA_LIBS $OPENMP_CXXFLAGS"

#############################################################
# Download + Build PVODE '98
#############################################################

CREATE_DIR(externalpackages)
CREATE_DIR(lib)
CREATE_DIR(include)

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
# Print configuration info
#############################################################

echo 
echo "Configuration summary"  | tee -a config-build.log
echo

if test "$FACETS_SOURCE" = ""
then
  echo "  FACETS support: no"   | tee -a config-build.log
else
  echo "  FACETS support: yes"  | tee -a config-build.log
fi

if test "$PETSC" = ""
then
  echo "  PETSc support: no" | tee -a config-build.log
else
  echo "  PETSc support: yes (version $PETSC_VERSION_MAJOR.$PETSC_VERSION_MINOR, release = $PETSC_VERSION_RELEASE)" | tee -a config-build.log
  echo "  PETSc has SUNDIALS support: $PETSC_HAS_SUNDIALS" | tee -a config-build.log
fi

if test "$IDALIBS" = ""
then
  echo "  IDA support: no" | tee -a config-build.log
else
  echo "  IDA support: yes" | tee -a config-build.log
fi

if test "$CVODELIBS" = ""
then
  echo "  CVODE support: no" | tee -a config-build.log
else
  echo "  CVODE support: yes" | tee -a config-build.log
fi

if test "$NCPATH" = ""
then
  echo "  NetCDF support: no" | tee -a config-build.log
else
  echo "  NetCDF support: yes" | tee -a config-build.log
fi

if test "$PNCPATH" = ""
then
  echo "  Parallel-NetCDF support: no" | tee -a config-build.log
else
  echo "  Parallel-NetCDF support: yes" | tee -a config-build.log
fi

if test "$PACT" = ""
then
  echo "  PDB support: no" | tee -a config-build.log
else
  echo "  PDB support: yes" | tee -a config-build.log
fi

if test "$HYPREPATH" = ""
then
  echo "  Hypre support: no" | tee -a config-build.log
else
  echo "  Hypre support: yes" | tee -a config-build.log
fi

if test "$MUMPS" = ""
then
  echo "  MUMPS support: no" | tee -a config-build.log
else
  echo "  MUMPS support: yes" | tee -a config-build.log
fi

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
