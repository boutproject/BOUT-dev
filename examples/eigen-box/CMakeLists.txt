cmake_minimum_required(VERSION 3.13)

project(eigen-box LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(eigen-box
  SOURCES eigen-box.cxx
  EXTRA_FILES eigenvals.py)
