
#include "domain.hxx"

#include <iostream>

using namespace std;

int main() {
  
  Domain d(10,20);
  Domain *d2 = d.splitX(3);
  Domain *d3 = d2->splitY(5);
  
  for(Domain::iterator it = d.begin(); it != d.end(); it++)
    cout << *it;
  
  return 0;
}
