
class BoutException;

#ifndef __BOUT_EXCEPTION_H__
#define __BOUT_EXCEPTION_H__

#include <exception>
#include <string>

using namespace std;

class BoutException : public exception {
public:
  BoutException(const char *, ...);
  virtual ~BoutException() throw();
  
  const char* what() const throw();

protected:
  string message;
};

#endif
