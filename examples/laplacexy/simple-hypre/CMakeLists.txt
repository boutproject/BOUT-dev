cmake_minimum_required(VERSION 3.13)

project(test_laplacexy LANGUAGES CXX C)
#find_package(bout++ REQUIRED)

if(NOT BOUT_HAS_HYPRE)
  message(FATAL_ERROR "This example requires BOUT++ to have been built with HYPRE")
endif()

set(test_laplacexy_hypre_sources
   test-laplacexy-hypre.cxx
)

set(test_laplacexy_hypre_depends bout++ PETSc HYPRE)

add_executable(test_laplacexy_hypre test-laplacexy-hypre.cxx)
target_link_libraries(test_laplacexy_hypre PRIVATE bout++::bout++)
target_compile_definitions(test_laplacexy_hypre PUBLIC "BOUT_HAS_HYPRE")
target_include_directories(test_laplacexy_hypre PUBLIC 
      ${HYPRE_DIR}/include
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
      ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
      )
if(BOUT_USE_CUDA)
   list(APPEND test_laplacexy_hypre_depends cuda)
   set_source_files_properties(${test_laplacexy_hypre_sources} PROPERTIES LANGUAGE CUDA )
   message(STATUS "test_laplacexy_hypre needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}") 
   target_compile_definitions(test_laplacexy_hypre PUBLIC "HYPRE_USING_CUDA;HYPRE_USING_UNIFIED_MEMORY")
   target_include_directories(test_laplacexy_hypre PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
endif ()

message(STATUS "test_laplacexy_hypre_depends: ${test_laplacexy_hypre_depends}")


  

