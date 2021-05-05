# Creates "bout++-time.cxx" in the build directory with the
# compilation date and time as variables

string(TIMESTAMP bout_date "%b %d %Y")
string(TIMESTAMP bout_time "%H:%M:%S")

set(bout_date_time_file
  "const char* boutcompiledate{\"${bout_date}\"}; const char* boutcompiletime{\"${bout_time}\"};")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/bout++-time.cxx" "${bout_date_time_file}")
