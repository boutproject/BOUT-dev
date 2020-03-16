cmake_minimum_required(VERSION 3.13)

project(conduction-snb LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(conduction-snb conduction-snb.cxx)
target_link_libraries(conduction-snb PRIVATE bout++::bout++)
