
class BoutException;

#ifndef __BOUT_EXCEPTION_H__
#define __BOUT_EXCEPTION_H__

#include <exception>
#include <string>

using std::string;

void BoutParallelThrowRhsFail(int status, const char* message);

class BoutException : public std::exception {
public:
  BoutException(const char *, ...);
  BoutException(const std::string);
  virtual ~BoutException();
  
  const char* what() const noexcept override;
  void Backtrace();
protected:
  char *buffer = nullptr;
  static const int BUFFER_LEN = 1024; // Length of char buffer for printing
  int buflen; // Length of char buffer for printing
  string message;
#ifdef BACKTRACE
  static const unsigned int TRACE_MAX = 128;
  void *trace[TRACE_MAX];
  char **messages;
  int trace_size;
  mutable std::string _tmp;
#endif
  std::string BacktraceGenerate() const;
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
