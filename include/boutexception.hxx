
class BoutException;

#ifndef __BOUT_EXCEPTION_H__
#define __BOUT_EXCEPTION_H__

#include <exception>
#include <string>

#include "string_utils.hxx"

using std::string;

void BoutParallelThrowRhsFail(int &status, const std::string &message);

class BoutException : public std::exception {
public:
  template <typename... Args>
  BoutException(const std::string &message, Args... args) :
    BoutException(string_format(message, args...)) {}
  BoutException(const std::string &message);
  virtual ~BoutException() throw() {}
  
  const char* what() const throw() { return message.c_str(); }
  void Backtrace();
protected:
  std::string message;
};

class BoutRhsFail : public BoutException {
public:
  template <typename... Args>
  BoutRhsFail(const std::string &message, Args... args) :
    BoutException(string_format(message, args...)) {}
  BoutRhsFail(const std::string &message) :
    BoutException(message) {}
};

class BoutIterationFail : public BoutException {
public:
  template <typename... Args>
  BoutIterationFail(const std::string &message, Args... args) :
    BoutException(string_format(message, args...)) {}
  BoutIterationFail(const std::string &message) :
    BoutException(message) {}
};

#endif
