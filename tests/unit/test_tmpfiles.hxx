#pragma once

#include <filesystem>
#include <fmt/core.h>
#include <string>

namespace bout::testing {

/// Return a temporary directory that definitely exists
inline auto test_directory() {
  auto test_dir = std::filesystem::temp_directory_path() / "bout_tests";
  create_directory(test_dir);
  return test_dir;
}

/// Create a uniquely named temporary file that is automatically cleaned up
class TempFile {
  static int current_count() {
    static int count{0};
    return count++;
  }

  std::filesystem::path filename{test_directory()
                                 / fmt::format("tempfile_{}", current_count())};

public:
  TempFile() = default;
  TempFile(const TempFile&) = delete;
  TempFile(TempFile&&) = delete;
  TempFile& operator=(const TempFile&) = delete;
  TempFile& operator=(TempFile&&) = delete;

  ~TempFile() { std::filesystem::remove(filename); }

  // Enable conversions to std::string / const char*
  operator std::string() const { return filename.string(); }
  auto c_str() const { return filename.c_str(); }
  auto string() const { return filename.string(); }
};

} // namespace bout::testing
