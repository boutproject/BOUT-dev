cmake_minimum_required(VERSION 3.13)

project(conducting-wall-mode LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(conducting-wall-mode cwm.cxx)
target_link_libraries(conducting-wall-mode PRIVATE bout++::bout++)
