# FindnetCDFCxx
# ----------
#
# Find the netCDF C++ API
#
# This module uses the ``ncxx4-config`` helper script as a hint for
# the location of the NetCDF C++ library. It should be in your PATH.
#
# This module will define the following variables:
#
# ::
#
#   netCDFCxx_FOUND - true if netCDFCxx was found
#   netCDFCxx_VERSION - netCDFCxx version in format Major.Minor.Release
#   netCDFCxx_INCLUDE_DIRS - Location of the netCDFCxx includes
#   netCDFCxx_LIBRARIES - Required libraries
#
# This module will also export the ``netCDF::netcdf-cxx4`` target.
#
# You can also set the following variables:
#
# ``netCDFCxx_ROOT``
#   Specify the path to the netCDF C++ installation to use
#
# ``netCDFCxx_DEBUG``
#   Set to TRUE to get extra debugging output

include(BOUT++functions)

find_package(netCDFCxx QUIET CONFIG)
if (netCDFCxx_FOUND)
  set(netCDFCxx_FOUND TRUE)
  if (NOT TARGET netCDF::netcdf-cxx4)
    bout_add_library_alias(netCDF::netcdf-cxx4 netcdf-cxx4)
  endif()
  return()
endif()

find_package(netCDF REQUIRED)

find_program(NCXX4_CONFIG "ncxx4-config"
  PATHS "${netCDFCxx_ROOT}"
  PATH_SUFFIXES bin
  DOC "Path to netCDF C++ config helper"
  NO_DEFAULT_PATH
  )

find_program(NCXX4_CONFIG "ncxx4-config"
  DOC "Path to netCDF C++ config helper"
  )

get_filename_component(NCXX4_CONFIG_TMP "${NCXX4_CONFIG}" DIRECTORY)
get_filename_component(NCXX4_CONFIG_LOCATION "${NCXX4_CONFIG_TMP}" DIRECTORY)
if (netCDFCxx_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NCXX4_CONFIG_LOCATION = ${NCXX4_CONFIG_LOCATION}")
endif()

bout_inspect_netcdf_config(NCXX4_HINTS_INCLUDE_DIR "${NCXX4_CONFIG}" "--includedir")
bout_inspect_netcdf_config(NCXX4_HINTS_PREFIX "${NCXX4_CONFIG}" "--prefix")

find_path(netCDF_CXX_INCLUDE_DIR
  NAMES netcdf
  DOC "netCDF C++ include directories"
  HINTS
    "${netCDF_C_INCLUDE_DIR}"
    "${NCXX4_HINTS_INCLUDE_DIR}"
    "${NCXX4_HINTS_PREFIX}"
    "${NCXX4_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "include"
  )
if (netCDFCxx_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " netCDF_CXX_INCLUDE_DIR = ${netCDF_CXX_INCLUDE_DIR}"
    " NCXX4_HINTS_INCLUDE_DIR = ${NCXX4_HINTS_INCLUDE_DIR}"
    " NCXX4_HINTS_PREFIX = ${NCXX4_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(netCDF_CXX_INCLUDE_DIR)

find_library(netCDF_CXX_LIBRARY
  NAMES netcdf_c++4 netcdf-cxx4
  DOC "netCDF C++ library"
  HINTS
    "${NCXX4_HINTS_INCLUDE_DIR}"
    "${NCXX4_HINTS_PREFIX}"
    "${NCXX4_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "lib" "lib64"
  )
if (netCDFCxx_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " netCDF_CXX_LIBRARY = ${netCDF_CXX_LIBRARY}"
    " NCXX4_HINTS_INCLUDE_DIR = ${NCXX4_HINTS_INCLUDE_DIR}"
    " NCXX4_HINTS_PREFIX = ${NCXX4_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(netCDF_CXX_LIBRARY)

bout_inspect_netcdf_config(_ncxx4_version "${NCXX4_CONFIG}" "--version")
if (_ncxx4_version)
  string(REGEX REPLACE "netCDF-cxx4 \([0-9]+\.[0-9]+\.[0-9]+\)" "\\1" netCDFCxx_VERSION "${_ncxx4_version}")
  message(STATUS "Found netCDFCxx version ${netCDFCxx_VERSION}")
else ()
  message(WARNING "Couldn't get NetCDF version")
endif()

unset(_ncxx4_version)
unset(_netcdf_version_major)
unset(_netcdf_version_minor)
unset(_netcdf_version_patch)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(netCDFCxx
  REQUIRED_VARS netCDF_CXX_LIBRARY netCDF_CXX_INCLUDE_DIR
  VERSION_VAR netCDFCxx_VERSION)

if (netCDFCxx_FOUND)
  set(netCDFCxx_INCLUDE_DIRS "${netCDF_CXX_INCLUDE_DIR}")
  set(netCDFCxx_LIBRARIES "${netCDF_CXX_LIBRARY}")

  if (NOT TARGET netCDF::netcdf-cxx4)
    add_library(netCDF::netcdf-cxx4 UNKNOWN IMPORTED)
    set_target_properties(netCDF::netcdf-cxx4 PROPERTIES
      IMPORTED_LINK_INTERFACE_LIBRARIES netCDF::netcdf
      IMPORTED_LOCATION "${netCDF_CXX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${netCDF_CXX_INCLUDE_DIR}")
  endif ()
endif ()
