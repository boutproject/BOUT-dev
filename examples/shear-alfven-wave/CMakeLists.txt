cmake_minimum_required(VERSION 3.13)

project(shear-alfven-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(shear-alfven-wave 2fluid.cxx)
target_link_libraries(shear-alfven-wave PRIVATE bout++::bout++)
