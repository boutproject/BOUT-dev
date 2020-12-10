cmake_minimum_required(VERSION 3.13)

project(constraints-laplace-dae LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(constraints-laplace-dae
  SOURCES laplace_dae.cxx
  EXTRA_FILES simple_xz.nc)
