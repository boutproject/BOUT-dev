cmake_minimum_required(VERSION 3.13)

project(drift-wave-constraint LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(drift-wave-constraint SOURCES test-drift.cxx)
