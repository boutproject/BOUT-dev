cmake_minimum_required(VERSION 3.13)

project(test_laplacexy_hypre LANGUAGES CXX C)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

if(NOT BOUT_HAS_HYPRE)
  message(FATAL_ERROR "This example requires BOUT++ to have been built with HYPRE")
endif()

bout_add_example(test_laplacexy_hypre
  SOURCES test-laplacexy-hypre.cxx
  EXTRA_FILES data/BOUT.inp)

if(BOUT_USE_CUDA)
  set_source_files_properties(test-laplacexy-hypre.cxx PROPERTIES LANGUAGE CUDA )
endif()
