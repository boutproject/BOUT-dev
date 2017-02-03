
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
  
  for(const auto& it : d) {
    cout << it;
  
  d.splitY(5);
  
  cout << "Split Y\n";
  for(const auto& it : d) {
    cout << it;
  
  cout << "\n\n";
  Domain t(8,5);
  cout << t;
  
  partition(&t, 5);
  
  cout << "\nPartitioned\n";
  for(const auto& it : t) {
    cout << it;
  
  return 0;
}
