# Creates "bout++-time.cxx" in the build directory with the
# compilation date and time as variables

set(bout_date_time_file
  "const char* boutcompiledate{__DATE__}; const char* boutcompiletime{__TIME__};")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/bout++-time.cxx" "${bout_date_time_file}")
