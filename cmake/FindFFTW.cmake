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

if (FFTW_INCLUDE_DIRS)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDE_DIRS)

find_path (FFTW_INCLUDE_DIRS fftw3.h)

find_library (FFTW_LIBRARIES NAMES fftw3)

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
