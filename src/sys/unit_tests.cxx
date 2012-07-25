#include <bout/sys/range.hxx>

#include <boutexception.hxx>
#include <iostream>
#include <assert.h>

using namespace std;

/// Skeleton implementation of BoutException
BoutException::~BoutException() throw() {}
BoutException::BoutException(const char* s, ...) {}
const char* BoutException::what() const throw() {std::cout << "ERROR\n";}

int main() {
  RangeIterator it(1,4, RangeIterator(6,9));
  
  int sum = 0;
  for(it.first(); !it.isDone(); it.next()) {
    cout << *it;
    sum += *it;
  }
  
  cout << " sum = " << sum << endl;
  assert(sum == 40);

  // A more standard syntax for iteration
  int prod = 1;
  for(it.first(); it != RangeIterator::end(); it++) {
    cout << *it;
    prod *= *it;
  }
  cout << " prod = " << prod << endl;
  assert(prod == 72576);
  
  // Check that an empty range is done
  RangeIterator nullit(5, 4);
  assert(nullit.isDone());
  
  // Check initialisation without call to first()
  sum = 0;
  for(RangeIterator it(2,5); !it.isDone(); it++) {
    sum += *it;
  }
  assert(sum == 14);

  // Range subtraction
  RangeIterator set(2, 10, RangeIterator(14,17));
  set -= RangeIterator(4,6);
  set -= RangeIterator(9,15);
  sum = 0;
  for(set.first(); set != RangeIterator::end(); set++) {
    cout << *set;
    sum += *set;
  }
  assert(sum == 53);
  
  cout << "\n";
  return 0;
}
