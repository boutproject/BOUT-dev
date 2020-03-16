cmake_minimum_required(VERSION 3.13)

project(conduction LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(conduction conduction.cxx)
target_link_libraries(conduction PRIVATE bout++::bout++)
