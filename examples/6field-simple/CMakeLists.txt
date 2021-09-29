cmake_minimum_required(VERSION 3.13)

project(6field-simple LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(elm_6f
  SOURCES elm_6f.cxx
  EXTRA_FILES cbm18_dens8.grid_nx68ny64.nc)
