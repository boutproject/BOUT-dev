find_path (PETSC_INCLUDE_DIRS petsc.h HINTS ${PETSC_DIR}/include)
find_library(PETSC_LIBRARIES NAMES petsc HINTS ${PETSC_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDE_DIRS)
message(STATUS "Third Party PETSC USED: ${PETSC_INCLUDE_DIRS}")
mark_as_advanced(PETSC_INCLUDE_DIRS PETSC_LIBRARIES)
