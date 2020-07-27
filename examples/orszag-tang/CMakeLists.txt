cmake_minimum_required(VERSION 3.13)

project(orszag-tang LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(orszag-tang
  SOURCES mhd.cxx
  EXTRA_FILES data/otv.grd.nc)
