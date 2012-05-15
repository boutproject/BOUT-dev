
#include "domain.hxx"
#include <boutexception.hxx>
#include <iostream>
#include <assert.h>

using namespace std;

/// Skeleton implementation of BoutException
BoutException::~BoutException() throw() {}
BoutException::BoutException(const char* s, ...) {}
const char* BoutException::what() const throw() {std::cout << "ERROR\n";}

int main() {
  
  Domain d(10,20);
  Domain *d2 = d.splitX(3);
  Domain *d3 = d2->splitY(5);
  
  int n = 0;
  for(Domain::iterator it = d.begin(); it != d.end(); it++) {
    cout << *it;
    n++;
  }
  assert(n == 3);
  
  return 0;
}
