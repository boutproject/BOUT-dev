cmake_minimum_required(VERSION 3.13)

project(hw3d LANGUAGES CXX C)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(hw3d
  SOURCES hw.cxx
  EXTRA_FILES README.md
              data/BOUT.inp)

if(BOUT_USE_CUDA)
  set_source_files_properties(hw.cxx PROPERTIES LANGUAGE CUDA)
endif()

