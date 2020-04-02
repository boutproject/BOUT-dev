find_path(SILO_INCLUDE_DIRS NAMES silo.h HINTS ${SILO_DIR}/include)
find_library(SILO_LIBRARIES NAMES silo siloh5 HINTS ${SILO_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SILO DEFAULT_MSG SILO_LIBRARIES SILO_INCLUDE_DIRS)

mark_as_advanced(SILO_INCLUDE_DIRS SILO_LIBRARIES)
