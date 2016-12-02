
class BoutException;

#ifndef __BOUT_EXCEPTION_H__
#define __BOUT_EXCEPTION_H__

#include <exception>
#include <string>

using namespace std;

void BoutParallelThrowRhsFail(int &status, const char* message);

class BoutException : public exception {
public:
  BoutException(const char *, ...);
  BoutException(const std::string);
  virtual ~BoutException() throw();
  
  const char* what() const throw();
  void Backtrace();
protected:
  char * buffer;
  static const int BUFFER_LEN = 1024; // Length of char buffer for printing
  int buflen; // Length of char buffer for printing
  string message;
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
