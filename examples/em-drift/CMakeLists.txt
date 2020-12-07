cmake_minimum_required(VERSION 3.13)

project(em-drift LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(em-drift
  SOURCES 2fluid.cxx
  EXTRA_FILES uedge.grd_beta.nc)
