cmake_minimum_required(VERSION 3.13)

project(shear-alfven-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(shear-alfven-wave
  SOURCES 2fluid.cxx
  EXTRA_FILES uedge.grd_Te_10.nc)
