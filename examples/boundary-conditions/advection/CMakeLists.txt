cmake_minimum_required(VERSION 3.13)

project(advection LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(advection
  SOURCES advection.cxx
  DATA_DIRS central-dirichlet 
            central-free 
            central-free-o3 
            upwind)
