
class BoutException;

#ifndef __BOUT_EXCEPTION_H__
#define __BOUT_EXCEPTION_H__

#include "bout/build_config.hxx"

#include <exception>
#include <string>
#include <utility>

#include "bout/deprecated.hxx"
#include "bout/format.hxx"

#include "fmt/core.h"

/// Throw BoutRhsFail with \p message if any one process has non-zero
/// \p status
void BoutParallelThrowRhsFail(int status, const char* message);

class BoutException : public std::exception {
public:
  BoutException(std::string msg) : message(std::move(msg)) { makeBacktrace(); }

  template <class S, class... Args>
  BoutException(const S& format, const Args&... args)
      : BoutException(fmt::format(format, args...)) {}

  ~BoutException() override;

  const char* what() const noexcept override {
    if (std::getenv("BOUT_SHOW_BACKTRACE") != nullptr) {
      getBacktrace();
      return (backtrace_message + "\n" + message).c_str();
    }

    return message.c_str();
  }
  void DEPRECATED(Backtrace()) {};

  /// Return the exception message along with the MsgStack and
  /// backtrace (if available)
  std::string getBacktrace() const;

  const std::string header{"====== Exception thrown ======\n"};

protected:
  std::string message;
#if BOUT_USE_BACKTRACE
  static constexpr unsigned int TRACE_MAX = 128;
  void* trace[TRACE_MAX];
  int trace_size;
  char** messages;
#endif
  std::string backtrace_message{};

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
