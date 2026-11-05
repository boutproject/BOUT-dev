include_guard(GLOBAL)

option(USE_ROCTX "Enable ROCTX profiling ranges" OFF)
option(USE_NVTX "Enable NVTX profiling ranges" OFF)

if(USE_ROCTX AND USE_NVTX)
  message(FATAL_ERROR "Enable only one of USE_ROCTX and USE_NVTX")
endif()

set(BOUT_PROFILING_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/../tools/profiling")

function(bout_enable_profiling profiling_target)
  # bout_add_model uses the plain link signature for executables.
  get_target_property(profiling_target_type ${profiling_target} TYPE)
  set(profiling_link_scope PRIVATE)
  if(profiling_target_type STREQUAL "EXECUTABLE")
    set(profiling_link_scope "")
  endif()
  if(USE_NVTX OR USE_ROCTX)
    target_include_directories(${profiling_target} PRIVATE "${BOUT_PROFILING_INCLUDE_DIR}")
  endif()
  if(USE_ROCTX)
    find_package(rocprofiler-sdk-roctx REQUIRED CONFIG)
    target_link_libraries(${profiling_target} ${profiling_link_scope} rocprofiler-sdk-roctx::rocprofiler-sdk-roctx)
    target_compile_definitions(${profiling_target} PRIVATE USE_ROCTX=1)
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
    target_include_directories(${profiling_target} PRIVATE ${NVTX3_INCLUDE_DIR})
    target_link_libraries(${profiling_target} ${profiling_link_scope} ${CMAKE_DL_LIBS})
    target_compile_definitions(${profiling_target} PRIVATE USE_NVTX=1)
  endif()
endfunction()
