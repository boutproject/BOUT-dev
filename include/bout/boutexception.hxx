#ifndef BOUT_EXCEPTION_H
#define BOUT_EXCEPTION_H

#include <exception>
#include <string>
#include <utility>

#include "fmt/base.h"
#include "fmt/core.h"

/// Throw BoutRhsFail with \p message if any one process has non-zero
/// \p status
void BoutParallelThrowRhsFail(int status, const char* message);

class BoutException : public std::exception {
public:
  BoutException(const BoutException&) = default;
  BoutException(BoutException&&) = delete;
  BoutException& operator=(const BoutException&) = default;
  BoutException& operator=(BoutException&&) = delete;
  BoutException(std::string msg);

  template <class... Args>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  BoutException(fmt::format_string<Args...> format, Args&&... args)
      : BoutException(fmt::vformat(format, fmt::make_format_args(args...))) {}

  ~BoutException() override;

  const char* what() const noexcept override { return message.c_str(); }

  /// Return the exception message along with the MsgStack and
  /// backtrace (if available)
  std::string getBacktrace() const;

  static void enableBacktrace() { show_backtrace = true; }
  static void disableBacktrace() { show_backtrace = false; }

private:
  std::string message;

  static bool show_backtrace;
};

class BoutRhsFail : public BoutException {
public:
  BoutRhsFail(std::string message) : BoutException(std::move(message)) {}
  template <class... Args>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  BoutRhsFail(fmt::format_string<Args...> format, Args&&... args)
      : BoutRhsFail(fmt::vformat(format, fmt::make_format_args(args...))) {}
};

class BoutIterationFail : public BoutException {
public:
  BoutIterationFail(std::string message) : BoutException(std::move(message)) {}
  template <class... Args>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  BoutIterationFail(fmt::format_string<Args...> format, Args&&... args)
      : BoutIterationFail(fmt::vformat(format, fmt::make_format_args(args...))) {}
};

#endif
