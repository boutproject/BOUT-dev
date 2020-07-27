cmake_minimum_required(VERSION 3.13)

project(advection-reaction LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(split_operator
  SOURCES split_operator.cxx
  EXTRA_FILES simple_xz.nc)
