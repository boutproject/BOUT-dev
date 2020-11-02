cmake_minimum_required(VERSION 3.13)


set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")

set(bd_sources
   blob2d.cxx
)

set(bd_depends bout++ PETSc)

if(HAVE_CUDA)
   list(APPEND bd_depends cuda nvToolsExt)
endif ()

if(HAVE_RAJA)
   list(APPEND bd_depends RAJA)
endif ()

message(STATUS "CUDA_TOOLKIT_ROOT_DIR : ${CUDA_TOOLKIT_ROOT_DIR}")
message(STATUS "bd_depends: ${bd_depends}")
blt_add_executable(
   NAME bd
   SOURCES ${bd_sources}
   INCLUDES 
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
     ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
     ${CUDA_TOOLKIT_ROOT_DIR}/include
   DEFINES 
      "BOUT_VERSION_STRING=\"${BOUT_FULL_VERSION}\""
      "BOUT_VERSION_DOUBLE=${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}"
   DEPENDS_ON
     ${bd_depends})


  #target_link_libraries(
  #   test_laplacexy_hypre PUBLIC 
  # ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse_static.a
  # ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcurand_static.a
  # ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libculibos.a
  # ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas_static.a
  # ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublasLt_static.a
  #)


