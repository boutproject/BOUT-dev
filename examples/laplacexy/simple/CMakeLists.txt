cmake_minimum_required(VERSION 3.13)

#project(test_laplacexy LANGUAGES CXX C)

set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")
#find_package(bout++ REQUIRED)
if(NOT BOUT_HAS_HYPRE)
  message(FATAL_ERROR "This example requires BOUT++ to have been built with HYPRE")
endif()
#add_executable(test_laplacexy test-laplacexy.cxx)
#target_link_libraries(test_laplacexy PRIVATE bout++::bout++)

set(test_laplacexy_sources
   test-laplacexy.cxx
)

set(test_laplacexy_depends bout++ PETSc HYPRE)

set(HYPRE_USING_CUDA_DEFINE "")
set(HYPRE_USING_CUBLAS_DEFINE "")
set(MEMORY_DEFINE "")
if(HAVE_CUDA)
   list(APPEND test_laplacexy_depends cuda)
   set(HYPRE_USING_CUDA_DEFINE "HYPRE_USING_CUDA")
   set(MEMORY_DEFINE "HYPRE_USING_UNIFIED_MEMORY")
   set(HYPRE_USING_CUBLAS_DEFINE "HYPRE_USING_CUBLAS")
endif ()

message(STATUS "CUDA_TOOLKIT_ROOT_DIR : ${CUDA_TOOLKIT_ROOT_DIR}")
message(STATUS "test_laplacexy_depends: ${test_laplacexy_depends}")
blt_add_executable(
   NAME test_laplacexy
   SOURCES ${test_laplacexy_sources}
   INCLUDES 
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
     ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
     ${CUDA_TOOLKIT_ROOT_DIR}/include
   DEFINES 
      "BOUT_VERSION_STRING=\"${BOUT_FULL_VERSION}\""
      "BOUT_VERSION_DOUBLE=${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}"
      BOUT_HAS_HYPRE
      ${HYPRE_USING_CUDA_DEFINE}
      ${HYPRE_USING_CUBLAS_DEFINE}
      ${MEMORY_DEFINE}
   DEPENDS_ON
     ${test_laplacexy_depends})


target_link_libraries(
   test_laplacexy PUBLIC 
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse_static.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcurand_static.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libculibos.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas_static.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublasLt_static.a
)

  
