
#include "domain.hxx"

#include <boutexception.hxx>

using std::list;

Domain::Domain(int NX, int NY) : nx(NX), ny(NY) {
  // Add boundaries
  Bndry b;
  
  b.side = xlow;
  b.start = 0;
  b.end = ny-1;
  b.neighbour = NULL;
  boundary.push_back(b);
  
  b.side = xhigh;
  boundary.push_back(b);
  
  b.side = ylow;
  b.end = nx-1;
  boundary.push_back(b);
  
  b.side = yhigh;
  boundary.push_back(b);
}

Domain::~Domain() {
  // Remove pointers to here from neighbours
  for(list<Bndry>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    if(it->neighbour != NULL)
      it->neighbour->removeNeighbour(this);
  }
}

void Domain::addNeighbour(Domain *d, 
                          BndrySide side, 
                          int first, int last, 
                          int shift) {
  
  // Find which boundary(s) this intersects with
  for(list<Bndry>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    if( (first <= it->end) && (last >= it->start) ) {
      // Overlaps
      
      if(it->neighbour != NULL) {
        // Already has a neighbour
        throw BoutException("Domain already has a neighbour in range %d -> %d", it->start, it->end);
      }
      
      if( (first == it->start) && (last == it->end) ) {
        // Just replace boundary
        it->neighbour = d;
      }else if(first == it->start) {
        // New domain is in first part
        Bndry b = *it;
        b.end = last; 
        b.neighbour = d;
        boundary.push_back(b);
        
        it->start = last+1;
      }else if(last == it->end) {
        // New domain is in second part
        Bndry b = *it;
        b.start = first;
        b.neighbour = d;
        boundary.push_back(b);
        
        it->end = first-1;
      }else if( (first > it->start) && (last < it->end) ){
        // Somewhere in the middle
        
        Bndry b = *it;
        b.start     = first;
        b.end       = last;
        b.neighbour = d;
        boundary.push_back(b);

        b.start     = last+1;
        b.end       = it->end;
        b.neighbour = NULL;
        boundary.push_back(b);
        
        it->end = first-1; // Becomes first of three
      }else {
        throw BoutException("Invalid range for new neighbour");
      }
      return;
    }
  }
  throw BoutException("Range out of bounds");
}

void Domain::removeNeighbour(Domain *d) {
  for(list<Bndry>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    if(it->neighbour == d) {
      removeNeighbour(it);
    }
  }
}

void Domain::removeNeighbour(Domain *d, BndrySide side) {
  for(list<Bndry>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    if((it->neighbour == d) && (it->side == side)) {
      removeNeighbour(it);
    }
  }
}

void Domain::removeNeighbour(list<Bndry>::iterator &it) {
  it->neighbour = NULL;
  // Merge with neighbours if also empty
}

Domain* Domain::splitX(int xind) {
  Domain *d = new Domain(xind, ny); // New domain on the left
  
  nx -= xind;
  
  for(list<Bndry>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    if(it->neighbour != NULL) {
      // Only interested in non-null boundaries
      
      if(it->side == xlow) { // Lower X boundaries
        // Move this to the new domain
        
        it->neighbour->removeNeighbour(this, reverse(it->side));
        
        it->neighbour->addNeighbour(d, reverse(it->side), it->start, it->end);
        d->addNeighbour(it->neighbour, it->side, it->start, it->end);
        
        removeNeighbour(it);
      }else if(it->side != xhigh) {
        // A Y boundary
        if(it->end < xind) {
          // Move to new domain
          it->neighbour->removeNeighbour(this,reverse(it->side));
          
          it->neighbour->addNeighbour(d, reverse(it->side), it->start, it->end);
          d->addNeighbour(it->neighbour, it->side, it->start, it->end);
          
          removeNeighbour(it);
          
        }else if((it->start <= xind) &&
                 (it->end >= xind)) {
          // Crosses cut
          
        }
      }
    }
  }
  
  // Add boundary between domains
  d->addNeighbour(this, xhigh, 0, ny-1);
  addNeighbour(d, xlow, 0, ny-1);
  
  return d;
}

Domain* Domain::splitY(int yind) {
  Domain *d = new Domain(nx, yind);
  
  ny -= yind;

  return d;
}

Domain::BndrySide Domain::reverse(const BndrySide &side) {
  switch(side) {
  case xlow:  return xhigh;
  case xhigh: return xlow;
  case ylow:  return yhigh;
  case yhigh: return ylow;
  }
  throw BoutException("Invalid BndrySide passed to reverse()");
}

#ifdef UNIT
// Unit test

#include <iostream>

/// Skeleton implementation of BoutException
BoutException::~BoutException() throw() {}
BoutException::BoutException(const char* s, ...) {}
const char* BoutException::what() const throw() {std::cout << "ERROR\n";}

int main() {
  
  Domain d(10,10);
  
  
  
  return 0;
}

#endif
