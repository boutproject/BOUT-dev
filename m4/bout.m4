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
  save_LIBS=$LIBS
  save_LDFLAGS=$LDFLAGS
  save_CPPFLAGS=$CPPFLAGS
  AC_MSG_CHECKING([for lib$1])
  AC_LANG_PUSH([C++])
  BACL_found=no

  # Try with no extra libraries first
  AS_IF([test ."$5" = .yes], [extra_prefix=""], [extra_prefix="$5"])
  LIBS="$save_LIBS $EXTRA_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    extern "C"
    char $2();
  ]], [[return $2();]])],
  [BACL_found=yes
   BOUT_MSG_DEBUG([found $1 without path or library flag])],
  [])
  LIBS=$save_LIBS

  # Now try with explicitly linking library
  AS_IF([test $BACL_found != yes], [
    LIBS="$save_LIBS $EXTRA_LIBS -l$1"
    AS_IF([test ."$5" = .yes], [extra_prefix=""], [extra_prefix="$5"])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[
      extern "C"
      char $2();
    ]], [[return $2();]])],
    [BACL_found=yes
     EXTRA_LIBS="$EXTRA_LIBS -l$1"
     BOUT_MSG_DEBUG([found $1 without path])],
    [])
  ])

  AS_IF([test $BACL_found != yes], [
    for search_prefix in $extra_prefix /usr /opt $HOME $HOME/local /usr/local ; do
      for path in $search_prefix $search_prefix/lib $search_prefix/lib64 $search_prefix/x86_64-linux-gnu
      do
        AS_IF([test -d $path], [
          LIBS="$save_LIBS $EXTRA_LIBS -l$1"
          LDFLAGS="$save_LDFLAGS -L$path"
          BOUT_MSG_DEBUG([try link $1 with $path])
          AC_LINK_IFELSE([AC_LANG_PROGRAM([[
            extern "C"
            char $2();
          ]], [[return $2();]])],
          [BACL_found=yes
           EXTRA_LIBS="$EXTRA_LIBS -L$path -l$1"
           BOUT_MSG_DEBUG([found $1 with $path])
           break],
          [])
        ])
      done
      AS_IF([test .$BACL_found = .yes],break;)
    done
   ])

   AS_IF([test $BACL_found = yes], [
     AC_MSG_RESULT(yes)
   ], [
     AC_MSG_RESULT(no)
   ])

   AS_IF([test $BACL_found = yes], [$3],[$4])

   LIBS=$save_LIBS
   LDFLAGS=$save_LDFLAGS
   CPPFLAGS=$save_CPPFLAGS
   AC_LANG_POP([C++])
])


dnl define the macro to check for header files
dnl first  argument is the name of the header file
dnl secound  arg is to be executed if found
dnl third  arg is to be executed if not found
dnl forth  arg is an additional path to check
AC_DEFUN([BOUT_ADDPATH_CHECK_HEADER],[
  AC_MSG_CHECKING([for $1])

  save_CPPFLAGS=$CPPFLAGS
  BACH_found=no

  AS_IF([test ."$4" != .yes], [extra_prefix="$4"], [extra_prefix=""])

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
    #include <$1>
  ])], [BACH_found=yes
        BOUT_MSG_DEBUG([found $1 without path])
        break])

  AS_IF([test $BACH_found != yes], [
    for search_prefix in $extra_prefix /usr /opt $HOME $HOME/local /usr/local
    do
      for path in $search_prefix $search_prefix/include
      do
        AS_IF([test -d $path], [
          CPPFLAGS="$save_CPPFLAGS -I$path"
          BOUT_MSG_DEBUG([try compile $1 with $path])
          AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([
              #include <$1>
              ], [])],
            [BACH_found=yes
             EXTRA_INCS="$EXTRA_INCS -I$path"
             BOUT_MSG_DEBUG([found $1 with $path])
             break])
        ])
      done
      AS_IF([test .$BACH_found = .yes], [break;])
    done
  ])

  AS_IF([test $BACH_found = yes], [
    AC_MSG_RESULT(yes)
  ], [
    AC_MSG_RESULT(no)
    CPPFLAGS=$save_CPPFLAGS
  ])
  AS_IF([test .$BACH_found = .yes], [$2], [$3])
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

dnl First argument is lower case module name
dnl Second argument is test program includes
dnl Third argument is test program main body
dnl Fourth argument is includes for test program to determine the integer type
dnl Fifth argument is main body for test program to determine the integer type 
AC_DEFUN([BOUT_FIND_SUNDIALS_MODULE],[

  dnl Slightly complicated as we have to deal with shell indirection
  AS_VAR_COPY([with_module], [with_$1])
  module_upper=m4_toupper($1)

  AC_MSG_NOTICE([Searching for SUNDIALS $module_upper library])
  AS_IF([test "$with_module" = "yes"], [
    # No path specified. Try using sundials-config
    AC_PATH_PROG([sundials_config], [sundials-config], [no], [$with_sundials$PATH_SEPARATOR$PATH])
    AS_IF([test "x$sundials_config" != xno], [
       AC_MSG_WARN(
         [Found sundials-config, this means your version of SUNDIALS is < 2.6, and probably won't work])
       sundials_module_includes=`$sundials_config -m $1 -t p -l c -s cppflags`
       sundials_module_libs=`$sundials_config -m $1 -t p -l c -s libs`
    ], [
       AC_MSG_WARN([No sundials-config available, no path given, will try compiling with $module_upper anyway])
       sundials_module_includes=""
       # Need to link to libsundials_ida, libsundials_cvode or libsundials_arkode
       sundials_module_libs="-lsundials_$1 -lsundials_nvecparallel $SUNDIALS_EXTRA_LIBS"
    ])
    AC_LANG_PUSH([C++])
    AC_MSG_CHECKING([if we can compile with SUNDIALS $module_upper])
    save_LIBS=$LIBS
    save_CXXFLAGS=$CXXFLAGS
    LIBS="$save_LIBS $sundials_module_libs"
    CXXFLAGS="$save_CXXFLAGS $sundials_module_includes"
    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM([
$2
         ], [$3])],
      [sundials_config_worked=yes],
      [sundials_config_worked=no])
    AC_MSG_RESULT([$sundials_config_worked])
    AS_IF([test $sundials_config_worked = yes], [
      AC_MSG_NOTICE([Using SUNDIALS $module_upper solver])
    ], [
      AC_MSG_FAILURE([Could not compile SUNDIALS $module_upper program, check your SUNDIALS version])
    ])
    LIBS=$save_LIBS
    CXXFLAGS="$save_CXXFLAGS"
    AC_LANG_POP([C++])
  ], [
    # Specified with path
    AC_MSG_NOTICE([Checking for $module_upper header files])

    # Check whether user supplied path to $module_upper install dir...
    AC_CHECK_FILES([$with_module/include/$1/$1.h
                    $with_module/include/$1/$1_spgmr.h
                    $with_module/include/$1/$1_bbdpre.h
                    $with_module/include/nvector/nvector_parallel.h
                    $with_module/include/sundials/sundials_types.h],
      [sundials_module_includes_found=yes
       sundials_module_includes_path=$with_module/include],
      [sundials_module_includes_found=no])
    AS_IF([test $sundials_module_includes_found = no], [
      # ...or path to $module_upper lib dir
      AC_CHECK_FILES([$with_module/../include/$1/$1.h
                      $with_module/../include/$1/$1_spgmr.h
                      $with_module/../include/$1/$1_bbdpre.h
                      $with_module/../include/nvector/nvector_parallel.h
                      $with_module/../include/sundials/sundials_types.h],
        [sundials_module_includes_found=yes
         sundials_module_includes_path=$with_module/../include],
        [sundials_module_includes_found=no])
    ])

    AS_IF([test $sundials_module_includes_found = no],
      [AC_MSG_FAILURE([Missing one or more $module_upper headers])])

    AC_MSG_NOTICE([Found $module_upper include path: $sundials_module_includes_path])

    # We've now got the include directory and can specify what libraries we need
    sundials_module_includes="-I$sundials_module_includes_path"
    sundials_module_libs="-lsundials_$1 -lsundials_nvecparallel $SUNDIALS_EXTRA_LIBS"

    # Try compiling something simple with a few different common paths
    save_LIBS=$LIBS
    save_LDFLAGS=$LDFLAGS
    save_CPPFLAGS=$CPPFLAGS
    AC_LANG_PUSH([C++])
    for sundials_module_lib_path in "$with_module" "$with_module/lib" "$with_module/lib64"
    do
      AC_MSG_CHECKING([if SUNDIALS $module_upper library path is $sundials_module_lib_path])
      LIBS="$save_LIBS $sundials_module_libs"
      LDFLAGS="$save_LDFLAGS -L$sundials_module_lib_path"
      CPPFLAGS="$save_CPPFLAGS $sundials_module_includes"
      AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([
$2
           ], [$3])],
        [sundials_module_lib_path_found=yes],
        [sundials_module_lib_path_found=no])
      AC_MSG_RESULT([$sundials_module_lib_path_found])
      AS_IF([test "x$sundials_module_lib_path_found" = "xyes"], [break])
      LIBS=$save_LIBS
      LDFLAGS=$save_LDFLAGS
      CPPFLAGS="$save_CPPFLAGS"
    done
    AC_LANG_POP([C++])

    SUNDIALS_MODULE_LDFLAGS="-L$sundials_module_lib_path"
  ])

  AS_IF([test "x$sundials_module_lib_path_found" = "xno"],
    [AC_MSG_FAILURE([Cannot compile $module_upper program])])

  # Compile in the $module_upper solver
  AC_MSG_NOTICE([=> $module_upper solver enabled])
  EXTRA_LIBS="$EXTRA_LIBS $SUNDIALS_MODULE_LDFLAGS $sundials_module_libs"
  EXTRA_INCS="$EXTRA_INCS $sundials_module_includes"
  CXXFLAGS="$CXXFLAGS -DBOUT_HAS_$module_upper"

  dnl The following is slightly complicated, but basically we use
  dnl AS_TR_SH to construct a shell variable from the variable
  dnl module_upper. This causes some shell indirection though, so we
  dnl then use AS_VAR_SET to actually assign the value we want to it
  AS_VAR_SET([AS_TR_SH([BOUT_HAS_$module_upper])], [yes])
  AS_VAR_SET([AS_TR_SH([${module_upper}LIBS])], ["$SUNDIALS_MODULE_LDFLAGS $sundials_module_libs"])
  AS_VAR_SET([AS_TR_SH([${module_upper}INCS])], ["$sundials_module_includes"])

  # Now we have successfully found the library, we need to determine
  # whether $module_upper uses int or long. Try to compile a simple
  # program to check
  save_CXXFLAGS=$CXXFLAGS
  AC_MSG_NOTICE(["checking $module_upper types..."])
  AC_LANG_PUSH([C++])

  for sundials_int_type in int long; do
    AC_MSG_CHECKING([$sundials_int_type])
    eval sundials_type_name=AS_TR_SH([${module_upper}INT])
    CXXFLAGS="$CXXFLAGS $sundials_module_includes -D$sundials_type_name=$sundials_int_type"
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
      AC_MSG_FAILURE([*** Cannot compile $module_upper with either long or int])
      ])
  AC_LANG_POP([C++])

  # We can now add that macro definition to the compilation flags
  CXXFLAGS="$save_CXXFLAGS -D$sundials_type_name=$sundials_int_type"
])
