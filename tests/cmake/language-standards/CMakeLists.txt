cmake_minimum_required(VERSION 3.17)
project(BoutLanguageStandards LANGUAGES CXX)
include("${BOUT_SOURCE_DIR}/cmake/BOUT++functions.cmake")

function(expect_property target property expected)
  get_target_property(actual ${target} ${property})
  if(NOT actual STREQUAL expected)
    message(
      FATAL_ERROR "${target}.${property}: expected ${expected}, got ${actual}"
    )
  endif()
endfunction()

# GPU properties can be checked without a GPU compiler. This is a configuration
# regression test, not a CUDA/HIP compilation or runtime test.
add_library(defaults INTERFACE)
bout_target_set_language_standards(defaults)
foreach(language IN ITEMS CXX CUDA HIP)
  expect_property(defaults ${language}_STANDARD 20)
  expect_property(defaults ${language}_STANDARD_REQUIRED ON)
  if(DEFINED CMAKE_${language}_STANDARD)
    message(FATAL_ERROR "Helper changed ${language} directory default")
  endif()
  set(CMAKE_${language}_STANDARD 23)
endforeach()

add_library(newer INTERFACE)
# Model an explicitly chosen target standard independently of CMake defaults.
foreach(language IN ITEMS CXX CUDA HIP)
  set_property(TARGET newer PROPERTY ${language}_STANDARD 23)
endforeach()
bout_target_set_language_standards(newer)
foreach(language IN ITEMS CXX CUDA HIP)
  expect_property(newer ${language}_STANDARD 23)
  if(NOT CMAKE_${language}_STANDARD EQUAL 23)
    message(FATAL_ERROR "Helper changed caller's ${language} default")
  endif()
endforeach()

foreach(old_standard IN ITEMS 98 17)
  add_library(older_${old_standard} INTERFACE)
  foreach(language IN ITEMS CXX CUDA HIP)
    set_property(
      TARGET older_${old_standard} PROPERTY ${language}_STANDARD
                                            ${old_standard}
    )
  endforeach()
  bout_target_set_language_standards(older_${old_standard})
  foreach(language IN ITEMS CXX CUDA HIP)
    expect_property(older_${old_standard} ${language}_STANDARD 20)
  endforeach()
endforeach()

add_library(unrelated INTERFACE)
foreach(language IN ITEMS CXX CUDA HIP)
  set_property(TARGET unrelated PROPERTY ${language}_STANDARD 17)
  expect_property(unrelated ${language}_STANDARD 17)
endforeach()

set(BOUT_HAS_HIP ON)
bout_target_enable_hip(newer)
expect_property(newer HIP_STANDARD 23)
expect_property(newer HIP_EXTENSIONS OFF)
expect_property(newer LINKER_LANGUAGE HIP)
bout_target_enable_hip(defaults)
expect_property(defaults HIP_STANDARD 20)

# Exercise SetupCompilers' HIP branch without requiring a ROCm installation.
# Only language detection is stubbed; the production configuration runs.
function(enable_language language)
  if(NOT language STREQUAL "HIP")
    message(FATAL_ERROR "Unexpected language detection: ${language}")
  endif()
  set(hip_detection_requested
      ON
      PARENT_SCOPE
  )
endfunction()
set(BOUT_ENABLE_MPI OFF)
set(BOUT_ENABLE_HIP ON)
set(BOUT_ENABLE_CUDA OFF)
set(BOUT_ENABLE_RAJA ON)
set(BOUT_ENABLE_UMPIRE ON)
set(CMAKE_HIP_STANDARD_REQUIRED OFF)
set(CMAKE_HIP_EXTENSIONS ON)
include("${BOUT_SOURCE_DIR}/cmake/SetupCompilers.cmake")
if(NOT hip_detection_requested
   OR NOT CMAKE_HIP_STANDARD EQUAL 23
   OR CMAKE_HIP_STANDARD_REQUIRED
   OR NOT CMAKE_HIP_EXTENSIONS
)
  message(FATAL_ERROR "HIP setup changed caller's directory defaults")
endif()

# Do not raise the project's supported compiler minimum just for this probe.
set(BOUT_HAS_HIP OFF)
if("cxx_std_23" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  set(probe_standard 23)
else()
  set(probe_standard 20)
endif()
set(CMAKE_CXX_STANDARD ${probe_standard})
add_executable(cpu_probe probe.cxx)
bout_target_enable_hip(cpu_probe)
expect_property(cpu_probe CXX_STANDARD ${probe_standard})
if(probe_standard EQUAL 23)
  target_compile_definitions(cpu_probe PRIVATE CHECK_CXX23)
endif()
message(STATUS "BOUT language standards regression checks passed")
