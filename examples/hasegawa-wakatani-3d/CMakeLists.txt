cmake_minimum_required(VERSION 3.13)

project(hw3d LANGUAGES CXX C)

if(NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

option(USE_ROCTX "Enable ROCTX profiling ranges in hw3d" OFF)
option(USE_NVTX "Enable NVTX profiling ranges in hw3d" OFF)

if(USE_ROCTX AND USE_NVTX)
  message(FATAL_ERROR "Enable only one of USE_ROCTX and USE_NVTX")
endif()

set(hw3d_source hw.cxx)
if(USE_ROCTX)
  set(hw3d_source hw-roctx.cxx)
elseif(USE_NVTX)
  set(hw3d_source hw-nvtx.cxx)
endif()

bout_add_example(hw3d SOURCES ${hw3d_source})

if(USE_NVTX OR USE_ROCTX)
  target_include_directories(
    hw3d PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../../tools/profiling"
  )
endif()

if(BOUT_HAS_CUDA)
  set_source_files_properties(${hw3d_source} PROPERTIES LANGUAGE CUDA)
endif()

if(USE_ROCTX)
  find_package(rocprofiler-sdk-roctx REQUIRED CONFIG)
  target_link_libraries(hw3d rocprofiler-sdk-roctx::rocprofiler-sdk-roctx)
  target_compile_definitions(hw3d PRIVATE USE_ROCTX=1)
elseif(USE_NVTX)
  find_path(
    NVTX3_INCLUDE_DIR nvtx3/nvToolsExt.h
    HINTS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
          ${CUDA_INCLUDE_DIRS}
          ${CUDAToolkit_INCLUDE_DIRS}
          ENV
          CUDA_HOME
          ENV
          CUDA_PATH
    PATH_SUFFIXES include
  )
  if(NOT NVTX3_INCLUDE_DIR)
    message(
      FATAL_ERROR
        "NVTX3 headers not found. Set NVTX3_INCLUDE_DIR to the directory containing nvtx3/nvToolsExt.h"
    )
  endif()
  target_include_directories(hw3d PRIVATE ${NVTX3_INCLUDE_DIR})
  target_link_libraries(hw3d ${CMAKE_DL_LIBS})
  target_compile_definitions(hw3d PRIVATE USE_NVTX=1)
endif()
