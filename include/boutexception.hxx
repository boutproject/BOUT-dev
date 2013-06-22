
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
  virtual ~BoutException() throw();
  
  const char* what() const throw();

protected:
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
