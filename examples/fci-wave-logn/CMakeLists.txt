cmake_minimum_required(VERSION 3.13)

project(fci-wave-logn LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(fci-wave-logn
  SOURCES fci-wave.cxx
  DATA_DIRS boundary div-integrate expanded)
