include(CheckCXXCompilerFlag)

# Add warning FLAGS to TARGET if the compiler supports them
function(target_enable_cxx_warning_if_supported TARGET)
  set(multiValueArgs FLAGS)
  cmake_parse_arguments(TARGET_ENABLE_WARNING "" "" "${multiValueArgs}" ${ARGN})

  foreach (WARNING_FLAG IN LISTS TARGET_ENABLE_WARNING_FLAGS)
    string(REPLACE "-" "_" WARNING_FLAG_STRIPPED ${WARNING_FLAG})

    # Note that gcc ignores unknown flags of the form "-Wno-warning"
    # for backwards compatibility. Therefore we need to add the
    # positive form as an additional flag which it will choke on (if
    # it doesn't exist). See: https://gcc.gnu.org/wiki/FAQ#wnowarning
    string(FIND ${WARNING_FLAG} "Wno-" NEGATIVE_FLAG_${WARNING_FLAG_STRIPPED})
    if (NEGATIVE_FLAG_${WARNING_FLAG_STRIPPED} EQUAL -1)
      set(IS_NEGATIVE_FLAG_${WARNING_FLAG_STRIPPED} FALSE)
    else()
      set(IS_NEGATIVE_FLAG_${WARNING_FLAG_STRIPPED} TRUE)
    endif()

    if (IS_NEGATIVE_FLAG_${WARNING_FLAG_STRIPPED})
      set(ORIGINAL_FLAG ${WARNING_FLAG})
      string(REPLACE "no-" "" WARNING_FLAG ${WARNING_FLAG})
      message(STATUS "Found negative flag: ${ORIGINAL_FLAG}\n"
        "   replaced with ${WARNING_FLAG}")
    endif()

    check_cxx_compiler_flag(${WARNING_FLAG} HAS_FLAG_${WARNING_FLAG_STRIPPED})

    if (IS_NEGATIVE_FLAG_${WARNING_FLAG_STRIPPED})
      set(WARNING_FLAG ${ORIGINAL_FLAG})
    endif()

    if (HAS_FLAG_${WARNING_FLAG_STRIPPED})
      message(STATUS "Warning flag is supported by compiler: ${WARNING_FLAG}")      

      target_compile_options(${TARGET} PRIVATE ${WARNING_FLAG})
    else()
      message(STATUS "Warning flag not supported by compiler: ${WARNING_FLAG}")
    endif()
  endforeach()
endfunction()
