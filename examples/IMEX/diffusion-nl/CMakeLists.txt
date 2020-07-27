cmake_minimum_required(VERSION 3.13)

project(diffusion-nl LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(diffusion-nl SOURCES diffusion-nl.cxx)
