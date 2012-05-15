
#include "partition.hxx"

#include <utils.hxx> // for ROUND
#include <boutexception.hxx>

#include <list>
using std::list;

#include <algorithm>
using std::swap;

void partition(Domain* d, int n) {
  // Check inputs
  if( (d == NULL) || (n < 1) )
    throw new BoutException("Invalid inputs to partition: %u, %d\n", d, n);
  
  int alpha = d->xSize();
  int beta  = d->ySize();
  
  bool swapped = false;
  if(beta > alpha) {
    swapped = true;
    swap(alpha, beta);
  }
  // Now satisfies beta <= alpha
  
  
  
}

void partitionAll(Domain* d, int n) {
  // Get a list of domains to partition
  list<Domain*> domains;
  vector<int> size;
  
  int total = 0;
  int ndomains = 0;
  for(Domain::iterator it = d->begin(); it != d->end(); it++) {
    domains.push_back(&(*it));
    // Get size of domain
    int area = it->area();
    size.push_back(area);
    total += area;
    ndomains++;
  }
  
  if(n < ndomains)
    throw new BoutException("Cannot partition %d domains into %d pieces", ndomains, n);

  // Now partition each one
  for(list<Domain*>::iterator it = domains.begin(); it != domains.end(); it++) {
    int area = (*it)->area();
    int nsub = ROUND(n * (((double) area) / total));
    if(nsub < 1)
      nsub = 1;
    partition(*it, nsub);
    total -= area; // Remaining area
    n -= nsub; // Remaining number of processors
  }
}


