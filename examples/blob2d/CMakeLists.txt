cmake_minimum_required(VERSION 3.13)

project(blob2d LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(blob2d
  SOURCES blob2d.cxx
  DATA_DIRS delta_0.25 delta_1 delta_10 two_blobs)
