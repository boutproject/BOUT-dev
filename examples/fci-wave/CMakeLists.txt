cmake_minimum_required(VERSION 3.13)

project(fci-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(fci-wave fci-wave.cxx)
target_link_libraries(fci-wave PRIVATE bout++::bout++)
