# FindNetCDF
# ----------
#
# Find the NetCDF IO library
#
# This module uses the ``nc-config`` and ``ncxx4-config`` helper scripts
# as hints for the location of the NetCDF libraries. They should be in
# your PATH.
#
# This module will define the following variables:
#
# ::
#
#   NetCDF_FOUND - true if NetCDF was found
#   NetCDF_VERSION - NetCDF version in format Major.Minor.Release
#   NetCDF_INCLUDE_DIRS - Location of the NetCDF includes
#   NetCDF_LIBRARIES - Required libraries
#
# This module will also export ``NetCDF::NetCDF_C`` and
# ``NetCDF::NetCDF_CXX`` targets.
#
# You can also set the following variables:
#
# ``NetCDF_ROOT``
#   Specify the path to the NetCDF installation to use
#
# ``NetCDF_DEBUG``
#   Set to TRUE to get extra debugging output


# Taken from https://github.com/conan-io/conan/issues/2125#issuecomment-351176653
# This is needed so we can make a clone of the NetCDF C++ target which
# has the name "netcdf-cxx4" by default
function(add_cloned_imported_target dst src)
    add_library(${dst} INTERFACE IMPORTED)
    foreach(name INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS)
        get_property(value TARGET ${src} PROPERTY ${name} )
        set_property(TARGET ${dst} PROPERTY ${name} ${value})
    endforeach()
endfunction()

find_package(netCDFCxx QUIET)
if (netCDFCxx_FOUND)
  if(NOT TARGET NetCDF::NetCDF_CXX)
    add_cloned_imported_target(NetCDF::NetCDF_CXX netCDF::netcdf-cxx4)
  endif()
  set(NetCDF_FOUND TRUE)
  return()
endif()

# A function to call nx-config with an argument, and append the resulting path to a list
# Taken from https://github.com/LiamBindle/geos-chem/blob/feature/CMake/CMakeScripts/FindNetCDF.cmake
function(inspect_netcdf_config VAR NX_CONFIG ARG)
    execute_process(
        COMMAND ${NX_CONFIG} ${ARG}
        OUTPUT_VARIABLE NX_CONFIG_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(EXISTS "${NX_CONFIG_OUTPUT}")
        set(${VAR} ${NX_CONFIG_OUTPUT} PARENT_SCOPE)
    endif()
endfunction()

find_program(NC_CONFIG "nc-config"
  PATHS "${NetCDF_ROOT}"
  PATH_SUFFIXES bin
  DOC "Path to NetCDF C config helper"
  NO_DEFAULT_PATH
  )

find_program(NC_CONFIG "nc-config"
  DOC "Path to NetCDF C config helper"
  )

get_filename_component(NC_CONFIG_TMP "${NC_CONFIG}" DIRECTORY)
get_filename_component(NC_CONFIG_LOCATION "${NC_CONFIG_TMP}" DIRECTORY)
if (NetCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NC_CONFIG_LOCATION = ${NC_CONFIG_LOCATION}"
    " NetCDF_ROOT = ${NetCDF_ROOT}")
endif()

find_program(NCXX4_CONFIG "ncxx4-config"
  PATHS "${NetCDF_ROOT}"
  PATH_SUFFIXES bin
  DOC "Path to NetCDF C++ config helper"
  NO_DEFAULT_PATH
  )

find_program(NCXX4_CONFIG "ncxx4-config"
  DOC "Path to NetCDF C++ config helper"
  )

get_filename_component(NCXX4_CONFIG_TMP "${NCXX4_CONFIG}" DIRECTORY)
get_filename_component(NCXX4_CONFIG_LOCATION "${NCXX4_CONFIG_TMP}" DIRECTORY)
if (NetCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NCXX4_CONFIG_LOCATION = ${NCXX4_CONFIG_LOCATION}")
endif()

inspect_netcdf_config(NC_HINTS_INCLUDE_DIR "${NC_CONFIG}" "--includedir")
inspect_netcdf_config(NC_HINTS_PREFIX "${NC_CONFIG}" "--prefix")

find_path(NetCDF_C_INCLUDE_DIR
  NAMES netcdf.h
  DOC "NetCDF C include directories"
  HINTS
    "${NC_HINTS_INCLUDE_DIR}"
    "${NC_HINTS_PREFIX}"
    "${NC_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "include"
  )
if (NetCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NetCDF_C_INCLUDE_DIR = ${NetCDF_C_INCLUDE_DIR}"
    " NC_HINTS_INCLUDE_DIR = ${NC_HINTS_INCLUDE_DIR}"
    " NC_HINTS_PREFIX = ${NC_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(NetCDF_C_INCLUDE_DIR)

find_library(NetCDF_C_LIBRARY
  NAMES netcdf
  DOC "NetCDF C library"
  HINTS
    "${NC_HINTS_INCLUDE_DIR}"
    "${NC_HINTS_PREFIX}"
    "${NC_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "lib" "lib64"
 )
if (NetCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NetCDF_C_LIBRARY = ${NetCDF_C_LIBRARY}"
    " NC_HINTS_INCLUDE_DIR = ${NC_HINTS_INCLUDE_DIR}"
    " NC_HINTS_PREFIX = ${NC_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(NetCDF_C_LIBRARY)

inspect_netcdf_config(NCXX4_HINTS_INCLUDE_DIR "${NCXX4_CONFIG}" "--includedir")
inspect_netcdf_config(NCXX4_HINTS_PREFIX "${NCXX4_CONFIG}" "--prefix")

find_path(NetCDF_CXX_INCLUDE_DIR
  NAMES netcdf
  DOC "NetCDF C++ include directories"
  HINTS
    "${NetCDF_C_INCLUDE_DIR}"
    "${NCXX4_HINTS_INCLUDE_DIR}"
    "${NCXX4_HINTS_PREFIX}"
    "${NCXX4_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "include"
  )
if (NetCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NetCDF_CXX_INCLUDE_DIR = ${NetCDF_CXX_INCLUDE_DIR}"
    " NCXX4_HINTS_INCLUDE_DIR = ${NCXX4_HINTS_INCLUDE_DIR}"
    " NCXX4_HINTS_PREFIX = ${NCXX4_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(NetCDF_CXX_INCLUDE_DIR)

find_library(NetCDF_CXX_LIBRARY
  NAMES netcdf_c++4 netcdf-cxx4
  DOC "NetCDF C++ library"
  HINTS
    "${NCXX4_HINTS_INCLUDE_DIR}"
    "${NCXX4_HINTS_PREFIX}"
    "${NCXX4_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "lib" "lib64"
  )
if (NetCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NetCDF_CXX_LIBRARY = ${NetCDF_CXX_LIBRARY}"
    " NCXX4_HINTS_INCLUDE_DIR = ${NCXX4_HINTS_INCLUDE_DIR}"
    " NCXX4_HINTS_PREFIX = ${NCXX4_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(NetCDF_CXX_LIBRARY)

if (NetCDF_INCLUDE_DIR)
  file(STRINGS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
    REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
  string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
  set(NetCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
  unset(_netcdf_version_major)
  unset(_netcdf_version_minor)
  unset(_netcdf_version_patch)
  unset(_netcdf_version_note)
  unset(_netcdf_version_lines)
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  REQUIRED_VARS NetCDF_C_LIBRARY NetCDF_C_INCLUDE_DIR NetCDF_CXX_LIBRARY NetCDF_CXX_INCLUDE_DIR
  VERSION_VAR NetCDF_VERSION)

if (NetCDF_FOUND)
  set(NetCDF_INCLUDE_DIRS "${NetCDF_CXX_INCLUDE_DIR}" "${NetCDF_C_INCLUDE_DIR}")
  set(NetCDF_LIBRARIES "${NetCDF_CXX_LIBRARY}" "${NetCDF_LIBRARY}")

  if (NOT TARGET NetCDF::NetCDF)
    add_library(NetCDF::NetCDF_C UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_C PROPERTIES
      IMPORTED_LOCATION "${NetCDF_C_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_C_INCLUDE_DIR}"
      )
    add_library(NetCDF::NetCDF_CXX UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_CXX PROPERTIES
      IMPORTED_LINK_INTERFACE_LIBRARIES NetCDF::NetCDF_C
      IMPORTED_LOCATION "${NetCDF_CXX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_CXX_INCLUDE_DIR}")
  endif ()
endif ()
