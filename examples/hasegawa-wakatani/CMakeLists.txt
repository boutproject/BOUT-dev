cmake_minimum_required(VERSION 3.13)

project(hasegama-wakatani LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(hasegama-wakatani hw.cxx)
target_link_libraries(hasegama-wakatani PRIVATE bout++::bout++)
