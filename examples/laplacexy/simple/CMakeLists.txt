cmake_minimum_required(VERSION 3.13)

project(laplacexy-simple LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(laplacexy-simple SOURCES test-laplacexy.cxx)
