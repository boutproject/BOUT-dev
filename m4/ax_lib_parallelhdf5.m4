# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_HDF5([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of HDF5 library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   former only looks for serial HDF5 installations via h5cc. The latter
#   only looks for parallel HDF5 installations via h5pcc. If the optional
#   argument is omitted, serial installations will be preferred over
#   parallel ones.
#
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the HDF5 library.
#     yes  - do check for HDF5 library in standard locations.
#     path - complete path to the HDF5 helper script h5cc or h5pcc.
#
#   If HDF5 is successfully found, this macro calls
#
#     AC_SUBST(HDF5_VERSION)
#     AC_SUBST(HDF5_CC)
#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)
#     AC_SUBST(HDF5_FC)
#     AC_SUBST(HDF5_FFLAGS)
#     AC_SUBST(HDF5_FLIBS)
#     AC_SUBST(HDF5_TYPE)
#     AC_DEFINE(HAVE_HDF5)
#
#   and sets with_hdf5="yes".  Additionally, the macro sets
#   with_hdf5_fortran="yes" if a matching Fortran wrapper script is found.
#   Note that Autconf's Fortran support is not used to perform this check.
#   H5CC and H5FC will contain the appropriate serial or parallel HDF5
#   wrapper script locations.
#
#   If HDF5 is disabled or not found, this macros sets with_hdf5="no" and
#   with_hdf5_fortran="no".
#
#   Your configuration script can test $with_hdf to take any further
#   actions. HDF5_{C,CPP,LD}FLAGS may be used when building with C or C++.
#   HDF5_F{FLAGS,LIBS} should be used when building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for HDF5 support
#        AX_LIB_HDF5()
#
#     2) dnl Check for serial HDF5 support
#        AX_LIB_HDF5([serial])
#
#     3) dnl Check for parallel HDF5 support
#        AX_LIB_HDF5([parallel])
#
#   One could test $with_hdf5 for the outcome or display it as follows
#
#     echo "HDF5 support:  $with_hdf5"
#
#   You could also for example, override the default CC in "configure.ac" to
#   enforce compilation with the compiler that HDF5 uses:
#
#     AX_LIB_HDF5([parallel])
#     if test "$with_hdf5" = "yes"; then
#             CC="$HDF5_CC"
#     else
#             AC_MSG_ERROR([Unable to find HDF5, we need parallel HDF5.])
#     fi
#
#   The HDF5_TYPE environment variable returns "parallel" or "serial",
#   depending on which type of library is found.
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 15

AC_DEFUN([AX_LIB_PARALLELHDF5], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Add a default --with-parallelhdf5 configuration option.
AC_ARG_WITH([parallelhdf5],
  AS_HELP_STRING(
    [--with-parallelhdf5=[yes/no/PATH]],
    [location of h5pcc for parallel HDF5 configuration]
  ),
  [if test "$withval" = "no"; then
     with_parallelhdf5="no"
   elif test "$withval" = "yes"; then
     with_parallelhdf5="yes"
   else
     with_parallelhdf5="yes"
     PARALLELH5CC="$withval"
   fi],
   [with_parallelhdf5="no"]
)

dnl Set defaults to blank
PARALLELHDF5_CC=""
PARALLELHDF5_VERSION=""
PARALLELHDF5_CFLAGS=""
PARALLELHDF5_CPPFLAGS=""
PARALLELHDF5_LDFLAGS=""
PARALLELHDF5_LIBS=""
PARALLELHDF5_FC=""
PARALLELHDF5_FFLAGS=""
PARALLELHDF5_FLIBS=""
PARALLELHDF5_TYPE=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_parallelhdf5" = "yes"; then
    if test -z "$PARALLELH5CC"; then
        dnl Check to see if H5CC is in the path.
        AC_PATH_PROGS(
            [PARALLELH5CC], [h5pcc]
            [])
    else
        AC_MSG_CHECKING([Using provided HDF5 C wrapper])
        AC_MSG_RESULT([$PARALLELH5CC])
    fi
    AC_MSG_CHECKING([for HDF5 type])
    AS_CASE([$PARALLELH5CC],
        [*h5pcc], [PARALLELHDF5_TYPE=parallel],
        [*h5cc], [PARALLELHDF5_TYPE=serial],
        [PARALLELHDF5_TYPE=neither])
    AC_MSG_RESULT([$PARALLELHDF5_TYPE])

    if test ! -f "$PARALLELH5CC" || test ! -x "$PARALLELH5CC"; then

        AC_MSG_CHECKING([if we can compile parallel HDF5 program without helper script])
        AC_LANG_PUSH([C++])
        AC_LINK_IFELSE(
          [AC_LANG_PROGRAM([
            #include <hdf5.h>
            ], [H5Fcreate(0, 0, 0, 0);])],
            [ac_cv_parallelhdf5_h=yes
             ac_cv_libparallelhdf5=yes],
            [ac_cv_parallelhdf5_h=no
             ac_cv_libparallelhdf5=no])
        AC_MSG_RESULT([$ac_cv_parallelhdf5_h])

        if test "$ac_cv_parallelhdf5_h" = "yes"; then
          AC_MSG_CHECKING([if found HDF5 is parallel])
          AC_EGREP_CPP([yes], [
            #include <H5pubconf.h>
            #ifdef H5_HAVE_PARALLEL
              yes
            #endif
          ], [ac_cv_parallelhdf5_h=yes], [ac_cv_parallelhdf5_h=no])
          AC_MSG_RESULT([$ac_cv_parallelhdf5_h])
        fi
        AC_LANG_POP([C++])

        if test "$ac_cv_parallelhdf5_h" = "no" ; then
          AC_MSG_FAILURE([
Unable to locate parallel HDF5 compilation helper script 'h5pcc'.
Please specify --with-parallelhdf5=<LOCATION> as the full path to h5pcc.
HDF5 support is being disabled (equivalent to --with-parallelhdf5=no).
])
          with_parallelhdf5="no"
        fi
    else
        AC_MSG_CHECKING([for HDF5 libraries])
        dnl Get the h5cc output
        PARALLELHDF5_SHOW=$(eval $PARALLELH5CC -show)

        dnl Get the actual compiler used
        PARALLELHDF5_CC=$(eval $PARALLELH5CC -show | $AWK '{print $[]1}')
        if test "$PARALLELHDF5_CC" = "ccache"; then
            PARALLELHDF5_CC=$(eval $PARALLELH5CC -show | $AWK '{print $[]2}')
        fi

        dnl h5cc provides both AM_ and non-AM_ options
        dnl depending on how it was compiled either one of
        dnl these are empty. Lets roll them both into one.

        dnl Look for "HDF5 Version: X.Y.Z"
        PARALLELHDF5_VERSION=$(eval $PARALLELH5CC -showconfig | $GREP 'HDF5 Version:' \
            | $AWK '{print $[]3}')

        dnl A ideal situation would be where everything we needed was
        dnl in the AM_* variables. However most systems are not like this
        dnl and seem to have the values in the non-AM variables.
        dnl
        dnl We try the following to find the flags:
        dnl (1) Look for "NAME:" tags
        dnl (2) Look for "H5_NAME:" tags
        dnl (3) Look for "AM_NAME:" tags
        dnl
        PARALLELHDF5_tmp_flags=$(eval $PARALLELH5CC -showconfig \
            | $GREP 'FLAGS\|Extra libraries:' \
            | $AWK -F: '{printf("%s "), $[]2}' )

        dnl Find the installation directory and append include/
        PARALLELHDF5_tmp_inst=$(eval $PARALLELH5CC -showconfig \
            | $GREP 'Installation point:' \
            | $AWK '{print $[]NF}' )

        dnl Add this to the CPPFLAGS
        PARALLELHDF5_CPPFLAGS="-I${PARALLELHDF5_tmp_inst}/include"

        dnl Now sort the flags out based upon their prefixes
        for arg in $PARALLELHDF5_SHOW $PARALLELHDF5_tmp_flags ; do
          case "$arg" in
            -I*) echo $PARALLELHDF5_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || PARALLELHDF5_CPPFLAGS="$arg $PARALLELHDF5_CPPFLAGS"
              ;;
            -L*) echo $PARALLELHDF5_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || PARALLELHDF5_LDFLAGS="$arg $PARALLELHDF5_LDFLAGS"
              ;;
            -l*) echo $PARALLELHDF5_LIBS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || PARALLELHDF5_LIBS="$arg $PARALLELHDF5_LIBS"
              ;;
          esac
        done

        PARALLELHDF5_LIBS="$PARALLELHDF5_LIBS -lhdf5"
        AC_MSG_RESULT([yes (version $[PARALLELHDF5_VERSION])])

        dnl See if we can compile
        AC_LANG_PUSH([C])
        ax_lib_parallelhdf5_save_CC=$CC
        ax_lib_parallelhdf5_save_CPPFLAGS=$CPPFLAGS
        ax_lib_parallelhdf5_save_LIBS=$LIBS
        ax_lib_parallelhdf5_save_LDFLAGS=$LDFLAGS
        CC=$PARALLELHDF5_CC
        CPPFLAGS=$PARALLELHDF5_CPPFLAGS
        LIBS=$PARALLELHDF5_LIBS
        LDFLAGS=$PARALLELHDF5_LDFLAGS
        AC_CHECK_HEADER([hdf5.h], [ac_cv_parallelhdf5_h=yes], [ac_cv_parallelhdf5_h=no])
        AC_CHECK_LIB([hdf5], [H5Fcreate], [ac_cv_libparallelhdf5=yes],
                     [ac_cv_libparallelhdf5=no])
        if test "$ac_cv_parallelhdf5_h" = "no" && test "$ac_cv_libparallelhdf5" = "no" ; then
          AC_MSG_FAILURE([Unable to compile HDF5 test program])
        fi
        dnl Look for HDF5's high level library
        AC_HAVE_LIBRARY([hdf5_hl], [PARALLELHDF5_LIBS="$PARALLELHDF5_LIBS -lhdf5_hl"], [], [])

        CC=$ax_lib_parallelhdf5_save_CC
        CPPFLAGS=$ax_lib_parallelhdf5_save_CPPFLAGS
        LIBS=$ax_lib_parallelhdf5_save_LIBS
        LDFLAGS=$ax_lib_parallelhdf5_save_LDFLAGS
        AC_LANG_POP([C])

	AC_SUBST([PARALLELHDF5_VERSION])
	AC_SUBST([PARALLELHDF5_CC])
	AC_SUBST([PARALLELHDF5_CFLAGS])
	AC_SUBST([PARALLELHDF5_CPPFLAGS])
	AC_SUBST([PARALLELHDF5_LDFLAGS])
	AC_SUBST([PARALLELHDF5_LIBS])
	AC_SUBST([PARALLELHDF5_TYPE])
	AC_DEFINE([HAVE_PARALLELHDF5], [1], [Defined if you have parallel HDF5 support])
    fi
fi
])
