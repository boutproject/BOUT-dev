find_path (PETSC_INCLUDE_DIRS petsc.h HINTS ${PETSC_DIR}/include)
find_library(PETSC_LIBRARIES NAMES petsc HINTS ${PETSC_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDE_DIRS)


if (NOT TARGET PETSc::PETSc)
  add_library(PETSc::PETSc UNKNOWN IMPORTED)
  list(GET PETSC_LIBRARIES 0 PETSC_LIBRARY)
  target_link_libraries(PETSc::PETSc INTERFACE "${PETSC_LIBRARIES}")
  set_target_properties(PETSc::PETSc PROPERTIES
    IMPORTED_LOCATION "${PETSC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIRS}"
    )
endif()
