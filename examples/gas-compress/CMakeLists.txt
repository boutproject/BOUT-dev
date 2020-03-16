cmake_minimum_required(VERSION 3.13)

project(gas-compress LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(gas-compress gas_compress.cxx gas_compress.hxx)
target_link_libraries(gas-compress PRIVATE bout++::bout++)
