#ifndef BOUT_EXCEPTION_H
#define BOUT_EXCEPTION_H

#include <exception>
#include <string>
#include <utility>

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

  template <class S, class... Args>
  BoutException(S&& format, Args&&... args)
      : BoutException(fmt::format(fmt::runtime(std::forward<S>(format)),
                                  std::forward<decltype(args)>(args)...)) {}

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
  template <class S, class... Args>
  BoutRhsFail(S&& format, Args&&... args)
      : BoutRhsFail(fmt::format(fmt::runtime(std::forward<S>(format)),
                                std::forward<decltype(args)>(args)...)) {}
};

class BoutIterationFail : public BoutException {
public:
  BoutIterationFail(std::string message) : BoutException(std::move(message)) {}
  template <class S, class... Args>
  BoutIterationFail(S&& format, Args&&... args)
      : BoutIterationFail(fmt::format(fmt::runtime(std::forward<S>(format)),
                                      std::forward<decltype(args)>(args)...)) {}
};

#endif
