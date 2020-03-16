cmake_minimum_required(VERSION 3.13)

project(em-drift LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(em-drift 2fluid.cxx)
target_link_libraries(em-drift PRIVATE bout++::bout++)
