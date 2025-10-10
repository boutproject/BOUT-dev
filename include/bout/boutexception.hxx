#ifndef BOUT_EXCEPTION_H
#define BOUT_EXCEPTION_H

#include "bout/build_defines.hxx"

#include <array>
#include <exception>
#include <string>
#include <utility>

#include "fmt/core.h"

/// Throw BoutRhsFail with \p message if any one process has non-zero
/// \p status
void BoutParallelThrowRhsFail(int status, const char* message);

class BoutException : public std::exception {
public:
  BoutException(std::string msg);

  template <class S, class... Args>
  BoutException(const S& format, const Args&... args)
      : BoutException(fmt::format(format, args...)) {}

  ~BoutException() override;

  const char* what() const noexcept override { return message.c_str(); }

  /// Return the exception message along with the MsgStack and
  /// backtrace (if available)
  std::string getBacktrace() const;

private:
  std::string message;
#if BOUT_USE_BACKTRACE
  static constexpr unsigned int TRACE_MAX = 128;
  std::array<void*, TRACE_MAX> trace{};
  int trace_size;
  char** messages;
#endif

  void makeBacktrace();
};

class BoutRhsFail : public BoutException {
public:
  BoutRhsFail(std::string message) : BoutException(std::move(message)) {}
  template <class S, class... Args>
  BoutRhsFail(const S& format, const Args&... args)
      : BoutRhsFail(fmt::format(format, args...)) {}
};

class BoutIterationFail : public BoutException {
public:
  BoutIterationFail(std::string message) : BoutException(std::move(message)) {}
  template <class S, class... Args>
  BoutIterationFail(const S& format, const Args&... args)
      : BoutIterationFail(fmt::format(format, args...)) {}
};

#endif
