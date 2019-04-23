if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)

find_library (FFTW_LIBRARIES NAMES fftw3)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)

if (FFTW_FOUND)
  if (NOT TARGET FFTW::FFTW)
    add_library(FFTW::FFTW UNKNOWN IMPORTED)
    set_target_properties(FFTW::FFTW PROPERTIES
      IMPORTED_LOCATION "${FFTW_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDES}"
      )
    # set_target_properties(FFTW::FFTW PROPERTIES
    #   IMPORTED_LOCATION "${FFTW_LIBRARIES}"
    #   INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
    #   )
  endif()
endif()
