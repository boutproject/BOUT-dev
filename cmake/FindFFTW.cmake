# FindFFTW
# ----------
#
# Find the Fastest Fourier Transform in the West FFT library
#
# This module uses the ``fftw-wisdom`` executable as a hint for the
# location of the FFTW library. It should be in your PATH.
#
# This module will define the following variables:
#
# ::
#
#   FFTW_FOUND - true if FFTW was found
#   FFTW_INCLUDE_DIRS - Location of the FFTW includes
#   FFTW_LIBRARIES - Required libraries
#
# This module will also export the ``FFTW::FFTW`` target.
#
# You can also set the following variables:
#
# ``FFTW_ROOT``
#   Specify the path to the FFTW installation to use
#
# ``FFTW_DEBUG``
#   Set to TRUE to get extra debugging output

if (FFTW_INCLUDE_DIRS)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDE_DIRS)

find_program(FFTW_WISDOM "fftw-wisdom"
  PATHS "${FFTW_ROOT}"
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH
  DOC "Path to fftw-wisdom executable"
  )

find_program(FFTW_WISDOM "fftw-wisdom"
  DOC "Path to fftw-wisdom executable"
  )
if (FFTW_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " FFTW_WISDOM = ${FFTW_WISDOM}"
    )
endif()

get_filename_component(FFTW_WISDOM_TMP "${FFTW_WISDOM}" DIRECTORY)
get_filename_component(FFTW_HINT_DIR "${FFTW_WISDOM_TMP}" DIRECTORY)

find_path(FFTW_INCLUDE_DIRS
  NAMES fftw3.h
  DOC "FFTW include directory"
  HINTS "${FFTW_HINT_DIR}"
  PATH_SUFFIXES "include"
  )
if (FFTW_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " FFTW_INCLUDE_DIRS = ${FFTW_INCLUDE_DIRS}"
    " FFTW_HINT_DIR = ${FFTW_HINT_DIR}"
    )
endif()

find_library (FFTW_LIBRARIES
  NAMES fftw3
  DOC "FFTW library location"
  HINTS "${FFTW_HINT_DIR}"
  PATH_SUFFIXES "lib" "lib64"
  )
if (FFTW_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " FFTW_LIBRARIES = ${FFTW_LIBRARIES}"
    " FFTW_HINT_DIR = ${FFTW_HINT_DIR}"
    )
endif()

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIRS)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDE_DIRS)

if (FFTW_FOUND AND NOT TARGET FFTW::FFTW)
  add_library(FFTW::FFTW UNKNOWN IMPORTED)
  set_target_properties(FFTW::FFTW PROPERTIES
    IMPORTED_LOCATION "${FFTW_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
    )
endif()
