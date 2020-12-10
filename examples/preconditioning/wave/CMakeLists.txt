cmake_minimum_required(VERSION 3.13)

project(preconditioning-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(preconditioning-wave SOURCES test_precon.cxx)
