
#include "partition.hxx"

#include <boutexception.hxx>

#include <list>
using std::list;

#include <vector>
using std::vector;

#include <algorithm>
using std::swap;

#include <cmath>

// Cast to double, to make code below less cluttered
double D(int x) {return (double) x; }

int roundi(double x) {
  return (x > 0.0) ? (int) (x + 0.5) : (int) (x - 0.5);
}

// Vertical equidissection into k m-strips, and j (m+1)-strips
void V(Domain* d, int m, int k, int j) {
  // Check inputs
  if( (d == NULL) || (m < 1) || (k < 0) || (j < 0) || ((k+j) < 1) )
    throw BoutException("Invalid inputs to partition V: %u, %d, %d, %d\n", d, m, k, j);
  
  int nx = d->xSize();
  int ny = d->ySize();
  
  // Work out sizes
  int nxm = roundi( D(nx) / (  D(k) + D(j)*D(m+1)/D(m) ) );
  if(nxm < 1)
    nxm = 1;
  
  int limit = (j>0) ? k : k-1;
  for(int x=0; x < limit; x++) {
    // Split into an m strip
    if(nxm == d->xSize())
      nxm--;
    Domain* strip = d->splitX(nxm);
    
    // Split strip into m pieces
    int nym = roundi( D(ny) / D(m) );
    if(nym < 1)
      nym = 1;
    for(int y=0;y<m-1;y++) {
      if(nym == strip->ySize())
        nym--;
      strip->splitY(nym);
    }
  }
  nxm = roundi( (D(m+1)/D(m)) * D(nx) / (  D(k) + D(j)*D(m+1)/D(m) ) );
  if(nxm < 1)
    nxm = 1;
  for(int x=0;x < j-1; x++) {
    // Split into an (m+1) strip
    if(nxm == d->xSize())
      nxm--;
    Domain* strip = d->splitX(nxm);
    
    // Split into (m+1) pieces
    int nym = roundi( D(ny) / D(m+1) );
    if(nym < 1)
      nym = 1;
    for(int y=0;y<m;y++) {
      if(nym == strip->ySize())
        nym--;
      strip->splitY(nym);
    }
  }
  
  // Split final strip
  int n = (j == 0) ? m : m+1;
  int nym = roundi( D(ny) / D(n) );
  if(nym < 1)
    nym = 1;
  for(int y=0;y<n-1;y++) {
    if(nym == d->ySize())
      nym--;
    d->splitY(nym);
  }
}

// Horizontal equidissection. Same as V, but interchange X and Y
void H(Domain* d, int m, int k, int j) {
  // Check inputs
  if( (d == NULL) || (m < 1) || (k < 0) || (j < 0) || ((k+j) < 1) )
    throw BoutException("Invalid inputs to partition H: %u, %d, %d, %d\n", d, m, k, j);
  
  int nx = d->xSize();
  int ny = d->ySize();
  
  // Work out sizes
  int nym = roundi( D(ny) / (  D(k) + D(j)*D(m+1)/D(m) ) );
  
  int limit = (j>0) ? k : k-1;
  for(int y=0; y < limit; y++) {
    // Split into an m strip
    if(nym == d->ySize())
      nym--;
    Domain* strip = d->splitY(nym);
    
    // Split strip into m pieces
    int nxm = roundi( D(nx) / D(m) );
    for(int x=0;x<m-1;x++) {
      if(nxm == strip->xSize())
        nxm--;
      strip->splitX(nxm);
    }
  }
  nym = roundi( (D(m+1)/D(m)) * D(ny) / (  D(k) + D(j)*D(m+1)/D(m) ) );
  for(int y=0;y < j-1; y++) {
    // Split into an (m+1) strip
    Domain* strip = d->splitY(nym);
    
    // Split into (m+1) pieces
    int nxm = roundi( D(nx) / D(m+1) );
    for(int x=0;x<m;x++) {
      if(nxm == strip->xSize())
        nxm--;
      strip->splitX(nxm);
    }
  }
  
  // Split final strip
  int n = (j == 0) ? m : m+1;
  int nxm = roundi( D(nx) / D(n) );
  for(int x=0;x<n-1;x++) {
    if(nxm == d->xSize())
      nxm--;
    d->splitX(nxm);
  }
}

void partition(Domain* d, int n) {
  // Check inputs
  if( (d == NULL) || (n < 1) )
    throw BoutException("Invalid inputs to partition: %u, %d\n", d, n);
  if(n == 1)
    return; // Nothing to do

  int alpha = d->xSize(); // Length of horizontal side
  int beta  = d->ySize(); // Length of vertical side
  
  bool swapped = false;
  if(beta > alpha) {
    swapped = true;
    swap(alpha, beta);
  }
  // Now satisfies beta <= alpha
  double val = D(beta) * n / alpha;
  
  // Find an integer s such that s(s-1) < val <= s(s+1)
  int s = (int) ( 0.5 + 0.5*sqrt(1. + 4.*val) );
  
  // Check value
  if( s*(s-1) > val )
    s--;
  if( s*(s+1) < val )
    s++;
  // Final check
  if( (s*(s-1) > val) || (s*(s+1) < val) )
    throw BoutException("Partition couldn't find s(s-1) < %e < s(s+1)", val);

  int r = n - floor( D(n) / D(s) ) * s;
  int t = floor( D(n) / D(s) ) - r;
  
  if( D(s+1)/D(t + r + 1) < D(beta)/D(alpha) ) {
    // V(s, t, r) is minimal equidissection
    if(swapped) {
      H(d, s, t, r);  // Swapped vertical <-> horizontal
    }else
      V(d, s, t, r);
  }else if( D(beta)/D(alpha) <= D(s-1)/D(t+r) ) {
    // V(s-1, s-r, t+2r+1-s) is a minimal equidissection
    if(swapped) {
      H(d, s-1, s-r, t+2*r+1-s);
    }else
      V(d, s-1, s-r, t+2*r+1-s);
  }else {
    // H(r+t, s-r, r) is a minimal equidissection
    if(swapped) {
      V(d, r+t, s-r, r);
    }else
      H(d, r+t, s-r, r);
  }
}

void partitionAll(Domain* d, int n) {
  if( (d == NULL) || (n < 1) )
    throw BoutException("Invalid inputs to partitionAll: %u, %d\n", d, n);
  
  // Get a list of domains to partition
  list<Domain*> domains;
  vector<int> size;
  
  int total = 0;
  int ndomains = 0;
  for(auto&& it : *d) {
  // for(Domain::iterator it = d->begin(); it != d->end(); it++) {
    domains.push_back(&it);
    // Get size of domain
    int area = it.area();
    size.push_back(area);
    total += area;
    ndomains++;
  }

  if(n < ndomains)
    throw BoutException("Cannot partition %d domains into %d pieces", ndomains, n);

  // Now partition each one
  for(const auto& it : domains) {
    int area = it->area();
    int nsub = roundi(n * (D(area) / total));
    if(nsub < 1)
      nsub = 1;
    partition(it, nsub);
    total -= area; // Remaining area
    n -= nsub; // Remaining number of processors
  }
}


