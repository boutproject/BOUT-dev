cmake_minimum_required(VERSION 3.13)

project(blob2d-laplacexz LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(blob2d-laplacexz SOURCES blob2d.cxx)
