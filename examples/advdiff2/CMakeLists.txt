cmake_minimum_required(VERSION 3.13)

project(advdiff2 LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(advdiff2 
  SOURCES header.hxx
          init.cxx
          rhs.cxx
  EXTRA_FILES slab.grd.nc)
