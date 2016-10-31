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
