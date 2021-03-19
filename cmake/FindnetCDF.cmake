# FindnetCDF
# ----------
#
# Find the netCDF IO library
#
# This module uses the ``nc-config`` helper script as a hint for the
# location of the netCDF libraries. It should be in your PATH.
#
# This module will define the following variables:
#
# ::
#
#   netCDF_FOUND - true if netCDF was found
#   netCDF_VERSION - netCDF version in format Major.Minor.Release
#   netCDF_INCLUDE_DIRS - Location of the netCDF includes
#   netCDF_LIBRARIES - Required libraries
#
# This module will also export the ``netCDF::netcdf`` target.
#
# You can also set the following variables:
#
# ``netCDF_ROOT``
#   Specify the path to the netCDF installation to use
#
# ``netCDF_DEBUG``
#   Set to TRUE to get extra debugging output

include(BOUT++functions)

find_package(netCDF QUIET CONFIG)
if (netCDF_FOUND)
  set(netCDF_FOUND TRUE)
  if (NOT TARGET netCDF::netcdf)
    bout_add_library_alias(netCDF::netcdf netcdf)
  endif()
  return()
endif()

find_program(NC_CONFIG "nc-config"
  PATHS "${netCDF_ROOT}"
  PATH_SUFFIXES bin
  DOC "Path to netCDF C config helper"
  NO_DEFAULT_PATH
  )

find_program(NC_CONFIG "nc-config"
  DOC "Path to netCDF C config helper"
  )

get_filename_component(NC_CONFIG_TMP "${NC_CONFIG}" DIRECTORY)
get_filename_component(NC_CONFIG_LOCATION "${NC_CONFIG_TMP}" DIRECTORY)
if (netCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NC_CONFIG_LOCATION = ${NC_CONFIG_LOCATION}"
    " netCDF_ROOT = ${netCDF_ROOT}")
endif()

bout_inspect_netcdf_config(NC_HINTS_INCLUDE_DIR "${NC_CONFIG}" "--includedir")
bout_inspect_netcdf_config(NC_HINTS_PREFIX "${NC_CONFIG}" "--prefix")

find_path(netCDF_C_INCLUDE_DIR
  NAMES netcdf.h
  DOC "netCDF C include directories"
  HINTS
    "${NC_HINTS_INCLUDE_DIR}"
    "${NC_HINTS_PREFIX}"
    "${NC_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "include"
  )
if (netCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " netCDF_C_INCLUDE_DIR = ${netCDF_C_INCLUDE_DIR}"
    " NC_HINTS_INCLUDE_DIR = ${NC_HINTS_INCLUDE_DIR}"
    " NC_HINTS_PREFIX = ${NC_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(netCDF_C_INCLUDE_DIR)

find_library(netCDF_C_LIBRARY
  NAMES netcdf
  DOC "netCDF C library"
  HINTS
    "${NC_HINTS_INCLUDE_DIR}"
    "${NC_HINTS_PREFIX}"
    "${NC_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "lib" "lib64"
 )
if (netCDF_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " netCDF_C_LIBRARY = ${netCDF_C_LIBRARY}"
    " NC_HINTS_INCLUDE_DIR = ${NC_HINTS_INCLUDE_DIR}"
    " NC_HINTS_PREFIX = ${NC_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(netCDF_C_LIBRARY)

if (netCDF_C_INCLUDE_DIR)
  file(STRINGS "${netCDF_C_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
    REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
  string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
  string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
  if (NOT _netcdf_version_note STREQUAL "")
    # Make development version compare higher than any patch level
    set(_netcdf_version_note ".99")
  endif()
  set(netCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
  unset(_netcdf_version_major)
  unset(_netcdf_version_minor)
  unset(_netcdf_version_patch)
  unset(_netcdf_version_note)
  unset(_netcdf_version_lines)
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(netCDF
  REQUIRED_VARS netCDF_C_LIBRARY netCDF_C_INCLUDE_DIR
  VERSION_VAR netCDF_VERSION)

if (netCDF_FOUND)
  set(netCDF_INCLUDE_DIR "${netCDF_C_INCLUDE_DIR}")
  set(netCDF_INCLUDE_DIRS "${netCDF_C_INCLUDE_DIR}")
  set(netCDF_LIBRARIES "${netCDF_C_LIBRARY}")

  if (NOT TARGET netCDF::netcdf)
    add_library(netCDF::netcdf UNKNOWN IMPORTED)
    set_target_properties(netCDF::netcdf PROPERTIES
      IMPORTED_LOCATION "${netCDF_C_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${netCDF_C_INCLUDE_DIR}"
      )
  endif ()
endif ()
