
#include "domain.hxx"
#include "partition.hxx"

#include <boutexception.hxx>
#include <iostream>
#include <assert.h>

using namespace std;

/// Skeleton implementation of BoutException
BoutException::~BoutException() throw() {}
BoutException::BoutException(const char* s, ...) {}
const char* BoutException::what() const throw() {std::cout << "ERROR\n";}

int main() {
  
  Domain d(20,10);
  d.splitX(5);
  
  cout << "Split X\n";
  
  for(Domain::iterator it = d.begin(); it != d.end(); it++)
    cout << *it;
  
  d.splitY(5);
  
  cout << "Split Y\n";
  for(Domain::iterator it = d.begin(); it != d.end(); it++)
    cout << *it;
  
  cout << "\n\n";
  Domain t(8,5);
  cout << t;
  
  partition(&t, 5);
  
  cout << "\nPartitioned\n";
  for(Domain::iterator it = t.begin(); it != t.end(); it++)
    cout << *it;
  
  return 0;
}
