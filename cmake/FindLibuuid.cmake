# FindLibuuid
# ----------
#
# Find a libuuid UUID generator library
#
# This module will define the following variables:
#
# ::
#
#   Libuuid_FOUND - true if Libuuid was found
#   Libuuid_INCLUDE_DIRS - Location of the Libuuid includes
#   Libuuid_LIBRARIES - Required libraries
#
# This module will also export the ``Libuuid::libuuid`` target.
#
# You can also set the following variables:
#
# ``Libuuid_ROOT``
#   Specify the path to the Libuuid installation to use
#
# ``Libuuid_DEBUG``
#   Set to TRUE to get extra debugging output

include (FindPackageHandleStandardArgs)

if (WIN32)
  find_package_handle_standard_args(Libuuid DEFAULT_MSG)
  return()
endif()

if (APPLE)
  find_library(CFLIB CoreFoundation)
  find_package_handle_standard_args(Libuuid DEFAULT_MSG CFLIB)
  mark_as_advanced(${CFLIB})

  if (Libuuid_FOUND AND NOT TARGET Libuuid::libuuid)
    add_library(Libuuid::libuuid UNKNOWN IMPORTED)
    set_target_properties(Libuuid::libuuid PROPERTIES
      IMPORTED_LOCATION ${CFLIB}
      )
  endif()
  return()
endif ()

find_path(Libuuid_INCLUDE_DIRS uuid/uuid.h)
find_library(Libuuid_LIBRARIES uuid)

set (SYMLINK_SYSTEM_UUID OFF
     CACHE
     BOOL
     "Make symlinks to the found libuuid (workaround for potential confilicts in LD_LIBRARY_PATH")
if (SYMLINK_SYSTEM_UUID)
  cmake_minimum_required(VERSION 3.14)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/uuid)
  file(CREATE_LINK
       ${Libuuid_INCLUDE_DIRS}/uuid/uuid.h
       ${CMAKE_CURRENT_BINARY_DIR}/include/uuid/uuid.h
       COPY_ON_ERROR SYMBOLIC)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
  file(CREATE_LINK
       ${Libuuid_LIBRARIES}
       ${CMAKE_CURRENT_BINARY_DIR}/lib/libuuid.so
       COPY_ON_ERROR SYMBOLIC)
  set(Libuuid_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/include/ CACHE PATH "" FORCE)
  set(Libuuid_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/lib/libuuid.so CACHE FILEPATH "" FORCE)
  set(Libuuid_ROOT ${CMAKE_CURRENT_BINARY_DIR} CACHE FILEPATH "" FORCE)
endif()

find_package_handle_standard_args(Libuuid DEFAULT_MSG Libuuid_LIBRARIES Libuuid_INCLUDE_DIRS)

mark_as_advanced(Libuuid_LIBRARIES Libuuid_INCLUDE_DIRS)

if (Libuuid_FOUND AND NOT TARGET Libuuid::libuuid)
  add_library(Libuuid::libuuid UNKNOWN IMPORTED)
  set_target_properties(Libuuid::libuuid PROPERTIES
    IMPORTED_LOCATION "${Libuuid_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${Libuuid_INCLUDE_DIRS}"
    )
endif()
