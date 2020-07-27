cmake_minimum_required(VERSION 3.13)

project(2Dturbulence_multigrid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(esel SOURCES esel.cxx)
