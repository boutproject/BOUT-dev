
class BoutException;

#ifndef __BOUT_EXCEPTION_H__
#define __BOUT_EXCEPTION_H__

#include <exception>
#include <string>

#include "bout/deprecated.hxx"

/// Throw BoutRhsFail with \p message if any one process has non-zero
/// \p status
void BoutParallelThrowRhsFail(int status, const char* message);

class BoutException : public std::exception {
public:
  BoutException(const char *, ...);
  BoutException(std::string msg) : message(std::move(msg)) {
    backtrace_message = makeBacktrace();
  }
  ~BoutException() override;

  const char* what() const noexcept override {
    return message.c_str();
  }
  void DEPRECATED(Backtrace()) {};

  /// Return the exception message along with the MsgStack and
  /// backtrace (if available)
  std::string getBacktrace() const;

  std::string makeBacktrace() const;

  const std::string header{"====== Exception thrown ======\n"};

protected:
  char *buffer = nullptr;
  static constexpr int BUFFER_LEN = 1024; // Length of char buffer for printing
  int buflen; // Length of char buffer for printing
  std::string message;
#ifdef BACKTRACE
  static constexpr unsigned int TRACE_MAX = 128;
#endif
  std::string backtrace_message{};
};

class BoutRhsFail : public BoutException {
public:
  BoutRhsFail(const char *, ...);
};

class BoutIterationFail : public BoutException {
public:
  BoutIterationFail(const char *, ...);
};

#endif
