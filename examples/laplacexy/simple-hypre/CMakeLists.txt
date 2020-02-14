cmake_minimum_required(VERSION 3.13)

project(test_laplacexy_hypre LANGUAGES CXX C)

find_package(bout++ REQUIRED)
if(NOT BOUT_HAS_HYPRE)
  message(FATAL_ERROR "This example requires BOUT++ to have been built with HYPRE")
endif()
add_executable(test_laplacexy_hypre test-laplacexy.cxx)
target_link_libraries(test_laplacexy_hypre PRIVATE bout++::bout++)
