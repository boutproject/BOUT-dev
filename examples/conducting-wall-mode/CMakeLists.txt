cmake_minimum_required(VERSION 3.13)

project(conducting-wall-mode LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(conducting-wall-mode SOURCES cwm.cxx)
