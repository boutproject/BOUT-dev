cmake_minimum_required(VERSION 3.13)

project(elm_pb LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(elm_pb_outerloop
  SOURCES elm_pb_outerloop.cxx
  EXTRA_FILES cbm18_dens8.grid_nx68ny64.nc
)

if(BOUT_USE_CUDA)
  set_source_files_properties(elm_pb_outerloop.cxx PROPERTIES LANGUAGE CUDA )
endif()

