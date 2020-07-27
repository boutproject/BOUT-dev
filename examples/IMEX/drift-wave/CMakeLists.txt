cmake_minimum_required(VERSION 3.13)

project(drift-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(test-drift SOURCES test-drift.cxx)
