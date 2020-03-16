cmake_minimum_required(VERSION 3.13)

project(6field-simple LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(elm_6f elm_6f.cxx)
target_link_libraries(elm_6f PRIVATE bout++::bout++)
