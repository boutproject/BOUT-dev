set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")


# UMPIRE
if (ENABLE_UMPIRE)
  find_package(UMPIRE REQUIRED)
  if(UMPIRE_FOUND)
     set (BOUT_HAS_UMPIRE ON)
  else()
     set(BOUT_HAS_UMPIRE OFF)
  endif()
endif ()

# RAJA
if (ENABLE_RAJA)
  find_package(RAJA REQUIRED)
  if (RAJA_FOUND)
     message(STATUS "RAJA_CONFIG:" ${RAJA_CONFIG}) 
     string(FIND ${RAJA_CONFIG} "raja" loc)
     math(EXPR value "${loc} + 5" OUTPUT_FORMAT DECIMAL)
     string(SUBSTRING ${RAJA_CONFIG} 0 ${value}  RAJA_PATH)
     message(STATUS "RAJA_PATH" ${RAJA_PATH})
     set(BOUT_HAS_RAJA ON)
  else()
     set(BOUT_HAS_RAJA OFF)   
  endif ()
endif ()


#HAVE_HYPRE
if (ENABLE_HYPRE OR HYPRE_DIR)
  enable_language(C)
  find_package(HYPRE REQUIRED)

  if(HYPRE_FOUND)
    set (BOUT_HAS_HYPRE ON)
    set (ENABLE_HYPRE ON) 

  endif ()
endif ()

# OpenMP
if (ENABLE_OPENMP)
  if (OPENMP_FOUND)
    set(HAVE_OPENMP True)
  endif ()
endif ()

#PETSc
if (ENABLE_PETSC OR PETSC_DIR)
  find_package(PETSc REQUIRED)

  if (PETSc_FOUND)
    set (BOUT_HAS_PETSC True)
    set (ENABLE_PETSC ON)

  endif ()
endif()
