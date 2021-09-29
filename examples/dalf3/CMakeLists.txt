cmake_minimum_required(VERSION 3.13)

project(dalf3 LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()


bout_add_example(dalf3
  SOURCES dalf3.cxx
  EXTRA_FILES cbm18_8_y064_x516_090309.nc)
