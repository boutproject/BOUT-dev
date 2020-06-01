find_path(HYPRE_INCLUDE_DIRS NAMES HYPRE.h HINTS ${HYPRE_DIR}/include)
find_library(HYPRE_LIBRARIES NAMES HYPRE HINTS ${HYPRE_DIR}/lib)

# Handle the QUIETLY and REQUIRED arguments and set HYPRE_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG HYPRE_LIBRARIES HYPRE_INCLUDE_DIRS)

mark_as_advanced(HYPRE_INCLUDE_DIRS HYPRE_LIBRARIES)
