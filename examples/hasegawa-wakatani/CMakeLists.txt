cmake_minimum_required(VERSION 3.13)

project(hasegawa-wakatani LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(hasegawa-wakatani SOURCES hw.cxx)

