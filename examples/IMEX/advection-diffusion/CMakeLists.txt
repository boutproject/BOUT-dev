cmake_minimum_required(VERSION 3.13)

project(advection-diffusion LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(imex SOURCES imex.cxx)
