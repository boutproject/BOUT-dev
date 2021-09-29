cmake_minimum_required(VERSION 3.13)

project(subsampling LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(subsampling
  SOURCES monitor.cxx
  EXTRA_FILES show.py)

