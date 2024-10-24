cmake_minimum_required(VERSION 3.13)

project(elm_pb_mri LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(elm_pb_mri
  SOURCES elm_pb_mri.cxx
  EXTRA_FILES data/BOUT.inp)
