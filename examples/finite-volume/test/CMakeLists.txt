cmake_minimum_required(VERSION 3.13)

project(finite-volume-test LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(finite-volume-test SOURCES finite_volume.cxx)
