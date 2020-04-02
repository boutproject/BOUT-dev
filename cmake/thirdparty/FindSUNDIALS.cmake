find_path(SUNDIALS_INCLUDE_DIRS NAMES sundials/sundials_config.h HINTS ${SUNDIALS_DIR}/include)

set(SUNDIALS_LIBS_LIST
  sundials_arkode 
  sundials_cvode
  sundials_cvodes
  sundials_farkode
  sundials_fcvode
  sundials_fida
  sundials_fkinsol
  sundials_fnvecparallel
  sundials_fnvecpthreads
  sundials_fnvecserial
  sundials_ida
  sundials_idas
  sundials_kinsol
  sundials_nvecparallel
  sundials_nvecpthreads
  sundials_nvecserial)

set(SUNDIALS_LIBRARIES)
foreach (LIB ${SUNDIALS_LIBS_LIST})
  find_library(SUNDIALS_LIB_${LIB}
    NAMES ${LIB}
    HINTS ${SUNDIALS_DIR}/lib)

  if (SUNDIALS_LIB_${LIB})
    set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARIES} ${SUNDIALS_LIB_${LIB}})
  endif ()
endforeach ()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SUNDIALS DEFAULT_MSG SUNDIALS_LIBRARIES SUNDIALS_INCLUDE_DIRS)

MARK_AS_ADVANCED(SUNDIALS_INCLUDE_DIRS SUNDIALS_LIBRARIES)
