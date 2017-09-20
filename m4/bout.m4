dnl -*- mode: autoconf -*-
dnl Copyright 2010-2016 B D Dudson, BOUT++ Team
dnl
dnl Contact Ben Dudson, bd512@york.ac.uk
dnl
dnl This file is part of BOUT++.
dnl
dnl BOUT++ is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU Lesser General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl BOUT++ is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public License
dnl along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
dnl
dnl ------------------------------------------------------------
dnl
dnl Additional macros for configure.ac
dnl Functions inspired by ESPResSo (espressomd.org)

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
    AC_LANG_PUSH([C++])
    LDFLAGS_save=$LDFLAGS
    LIBS="$EXTRA_LIBS -l$1"
    BACL_found=no
    AS_IF([test ."$5" = .yes], [extra_prefix=""],[extra_prefix="$5"])
    AC_TRY_LINK([extern "C"
    char $2();] ,[return $2();], [BACL_found=yes ; break;
              BOUT_MSG_DEBUG([found $1 without path])
        ],)
    if test $BACL_found != yes ; then
        for search_prefix in $extra_prefix /usr /opt $HOME $HOME/local /usr/local ; do
            for path in $search_prefix $search_prefix/lib $search_prefix/lib64 $search_prefix/x86_64-linux-gnu ; do
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
    AC_LANG_POP([C++])
])


dnl define the macro to check for header files
dnl first  argument is the name of the header file
dnl secound  arg is to be executed if found
dnl third  arg is to be executed if not found
dnl forth  arg is an additional path to check
AC_DEFUN([BOUT_ADDPATH_CHECK_HEADER],[
    AC_MSG_CHECKING([for $1])
    CPPFLAGS_save=$CPPFLAGS
    BACH_found=no
    AS_IF([test ."$4" != .yes], [extra_prefix="$4"],[extra_prefix=""])
    AC_TRY_COMPILE([#include <$1>] ,, [BACH_found=yes ; break;
            BOUT_MSG_DEBUG([found $1 without path])
        ],)
    if test $BACH_found != yes ; then
        for search_prefix in $extra_prefix /usr /opt $HOME $HOME/local /usr/local ; do
            for path in $search_prefix $search_prefix/include ; do
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

dnl Check the int type of packages in SUNDIALS
dnl First argument is name of the package (e.g. CVODE)
dnl Second argument is include path (e.g. $CVODEINCS)
dnl Third argument is macro definition (e.g. CVODEINT)
dnl Fourth argument is header and rhs function declaration
dnl Fifth argument is PrecInit function call
AC_DEFUN([BOUT_CHECK_SUNDIALS_TYPE],[
  # Try to compile a simple program to check whether $1 uses int or long
  save_CXXFLAGS=$CXXFLAGS
  AC_MSG_NOTICE(["checking $1 types..."])
  AC_LANG_PUSH([C++])

  for sundials_int_type in int long; do
    AC_MSG_CHECKING([$sundials_int_type])
    CXXFLAGS="$CXXFLAGS $2 -D$3=$sundials_int_type"
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([
$4
      ], [$5])],
    [sundials_int_type_found=yes],
    [sundials_int_type_found=no])
    AC_MSG_RESULT($sundials_int_type_found)
    if test "x$sundials_int_type_found" = "xyes"; then
      break;
    fi
    CXXFLAGS=$save_CXXFLAGS
  done

  AS_IF([test "x$sundials_int_type_found" = "xno"], [
      AC_MSG_FAILURE([*** Cannot compile $1 with either long or int])
      ])
  AC_LANG_POP([C++])
  CXXFLAGS="$save_CXXFLAGS -D$3=$sundials_int_type"
])

AC_DEFUN([BOUT_CHECK_PRETTYFUNCTION], [
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([does C++ compiler support __PRETTY_FUNCTION__])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[]],
                 [[const char* name = __PRETTY_FUNCTION__;]])],
    [AC_MSG_RESULT(yes)
     CXXFLAGS="$CXXFLAGS -DHAS_PRETTY_FUNCTION"],
    [AC_MSG_RESULT(no)])
  AC_LANG_POP([C++])
])
