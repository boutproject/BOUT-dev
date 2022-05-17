cmake_minimum_required(VERSION 3.13)

project(blob2d-outerloop LANGUAGES CXX C)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(blob2d-outerloop SOURCES blob2d.cxx)

if(BOUT_USE_CUDA)
  set_source_files_properties(blob2d.cxx PROPERTIES LANGUAGE CUDA)
endif()
