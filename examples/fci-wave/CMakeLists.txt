cmake_minimum_required(VERSION 3.13)

project(fci-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(fci-wave
  SOURCES fci-wave.cxx
  DATA_DIRS div div-integrate logn
  EXTRA_FILES compare-density.py)
