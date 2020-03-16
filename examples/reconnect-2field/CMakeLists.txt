cmake_minimum_required(VERSION 3.13)

project(reconnect-2field LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(reconnect-2field 2field.cxx)
target_link_libraries(reconnect-2field PRIVATE bout++::bout++)
