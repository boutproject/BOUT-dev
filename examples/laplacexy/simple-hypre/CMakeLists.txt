cmake_minimum_required(VERSION 3.13)

#project(test_laplacexy_hypre LANGUAGES CXX C)

#find_package(bout++ REQUIRED)
if(NOT BOUT_HAS_HYPRE)
  message(FATAL_ERROR "This example requires BOUT++ to have been built with HYPRE")
endif()
#add_executable(test_laplacexy_hypre test-laplacexy.cxx)
#target_link_libraries(test_laplacexy_hypre PRIVATE bout++::bout++)

set(test_laplacexy_hypre_sources
   test-laplacexy.cxx
)

set(test_laplacexy_hypre_depends bout++ PETSc HYPRE)

blt_add_executable(
   NAME test_laplacexy_hypre
   SOURCES ${test_laplacexy_hypre_sources}
   INCLUDES 
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
     ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
   DEFINES 
      "BOUT_VERSION_STRING=\"${BOUT_FULL_VERSION}\""
      "BOUT_VERSION_DOUBLE=${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}"
      BOUT_HAS_HYPRE
   DEPENDS_ON
     ${test_laplacexy_hypre_depends})


