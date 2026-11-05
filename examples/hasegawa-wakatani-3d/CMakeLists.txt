cmake_minimum_required(VERSION 3.13)

project(hw3d LANGUAGES CXX C)

if(NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/SetupProfiling.cmake")

set(hw3d_source hw.cxx)
if(USE_ROCTX)
  set(hw3d_source hw-roctx.cxx)
elseif(USE_NVTX)
  set(hw3d_source hw-nvtx.cxx)
endif()

bout_add_example(hw3d SOURCES ${hw3d_source})

bout_enable_profiling(hw3d)

if(BOUT_HAS_CUDA)
  set_source_files_properties(${hw3d_source} PROPERTIES LANGUAGE CUDA)
endif()
