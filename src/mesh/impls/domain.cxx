
#include "domain.hxx"

#include <boutexception.hxx>
#include <assert.h>

DomainIterator::DomainIterator() : dom(0) {}

DomainIterator::DomainIterator(Domain &d) : dom(&d) {
  stack.push_front(dom);
  visited_list.push_front(dom);
}

DomainIterator& DomainIterator::operator++() {
  // Depth-first search of domains, using a list of visited 
  // This algorithm is pretty slow, but should only be run during setup

  if(dom == NULL) // Already at the end
    return *this;

  while(!stack.empty()) {
    // Most recently visited first
    Domain* d = stack.front();
    
    // Go through all boundaries of this domain
    for(list<Domain::Bndry*>::iterator bit = d->boundary.begin(); bit != d->boundary.end(); bit++) {
      Domain::Bndry* b = *bit;
      Domain *neighbour = b->getNeighbour(d);
      if((neighbour != NULL) && !visited(neighbour)) {
        // If neighbour exists, and hasn't been visited
        stack.push_front(neighbour);
        visited_list.push_front(neighbour);
        dom = neighbour;
        return *this;
      }
    }
    // All boundaries visited
    stack.pop_front(); // Remove from stack, but not from visited_list
  }
  // No more domains
  dom = NULL;
  return *this;
}

Domain::Domain(int NX, int NY) : nx(NX), ny(NY) {}
Domain::Domain(int NX, int NY, const string &Name) : nx(NX), ny(NY), name(Name) {}

Domain::~Domain() {
  for(list<Bndry*>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    // Remove pointers to here from neighbours
    Domain::Bndry *b = *it;
    if(b->from == this)
      b->from = NULL;
    if(b->to == this)
      b->to = NULL;
    
    // If boundary doesn't connect to anything, delete
    if( (b->from == NULL) && (b->to == NULL) )
      delete b;
  }

  boundary.clear();
}

void Domain::addNeighbour(Domain *d, 
                          BndrySide side, 
                          int first, int last, 
                          int shift, int zshift) {
  // Create a new boundary
  Bndry *b = new Bndry(side, this, d, first, last, shift, zshift);
  
  // Add to both domains
  addBoundary(b);
  if((d != NULL) && (d != this))
    d->addBoundary(b);
}

#include <stdio.h>

void Domain::addBoundary(Bndry *b) {
  // Check that the range is within bounds
  BndrySide s = b->side;
  int len = nx;
  if( (s == xlow) || (s == xhigh) )
    len = ny;
  if( b->onSide(this, s) ) {
    if( (b->getMin(s) < 0) || (b->getMax(s) >= len) ) {
      fprintf(stderr, "%x: (%d, %d) %d  %d [%d, %d]",
              this, nx, ny, s, len, b->getMin(s), b->getMax(s));
      throw BoutException("Invalid range");
    }
  }
  s = reverse(s);
  if( b->onSide(this, s) ) {
    if( (b->getMin(s) < 0) || (b->getMax(s) >= len) ) {
      fprintf(stderr, "%x: (%d, %d) %d  %d [%d, %d]",
              this, nx, ny, s, len, b->getMin(s), b->getMax(s));
      throw BoutException("Invalid range");
    }
  }

  // Check that the range doesn't overlap an existing boundary
  for(list<Bndry*>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    if( ((*it)->side != b->side) && ((*it)->side != reverse(b->side)) )
      continue;
    
    BndrySide s = b->side;
    if((*it)->onSide(this, s) && b->onSide(this, s)) {
      // Same side as an existing boundary
      // Check for overlap
      if( !( ((*it)->getMax(s) < b->getMin(s)) || (b->getMax(s) < (*it)->getMin(s)) ) ) {
        fprintf(stderr, "%x: (%d, %d) %d [%d, %d] [%d, %d]",
                this, nx, ny, s, (*it)->getMin(s), (*it)->getMax(s), b->getMin(s), b->getMax(s));
        throw BoutException("Boundary ranges overlap");
      }
    }
    s = reverse(s);
    if((*it)->onSide(this, s) && b->onSide(this, s)) {
      // Same side as an existing boundary
      // Check for overlap
      if( !( ((*it)->getMax(s) < b->getMin(s)) || (b->getMax(s) < (*it)->getMin(s)) ) )
        throw BoutException("Boundary ranges overlap");
    }
  }
  
  // Add to list of boundaries
  boundary.push_back(b);
}

void Domain::removeBoundary(Bndry *b) {
  int s = boundary.size();
  boundary.remove(b);
  assert(boundary.size() == s-1);
}

void Domain::removeBoundary(list<Bndry*>::iterator &it) {
  it = boundary.erase(it);
}

#include <iostream>

Domain* Domain::splitX(int xind) {

  cout << "SplitX " << this << " (" << nx << "," << ny << ") " << xind << "\n" ;
  
  Domain *d = create(xind, ny); // New domain on the right
  nx -= xind;
  
  cout << "   " << this << " (" << nx << "," << ny << ")\n";
  cout << "   " << d << " (" << d->nx << "," << d->ny << ")\n";
  
  // Copy the list of boundaries so not iterating over a changing list
  list<Bndry*> oldboundary(boundary);
  
  for(list<Bndry*>::iterator it=oldboundary.begin(); it != oldboundary.end(); it++) {
    Bndry *b = *it;
    
    if(b->onSide(this, xhigh)) {
      // Move to new domain
      b->setNeighbour(xlow, d);
    }else if(b->onSide(this, ylow)) {
      int start = b->getMin(ylow);
      int end   = b->getMax(ylow);
      
      if(start >= nx) {
        // Move to new domain
        b->shiftInds(ylow, -nx);
        b->setNeighbour(yhigh, d);
      }else if(end >= nx) {
        // Crosses cut, so need to split into two
        
        b->end -= end - nx + 1;  // Make boundary smaller
        
        // Create a new boundary
        d->addNeighbour(b->getNeighbour(this), 
                        ylow,
                        0, end - nx,
                        b->getShift(ylow)+nx, b->zshift);
      }
    }else if(b->onSide(this, yhigh)) {
      int start = b->getMin(yhigh);
      int end   = b->getMax(yhigh);
      
      if(start >= nx) {
        // Move to new domain
        b->shiftInds(yhigh, -nx);
        b->setNeighbour(ylow, d);
      }else if(end >= nx) {
        // Crosses cut, so need to split into two
        b->end -= end - nx + 1;      // Make boundary smaller
        
        // Create a new boundary
        d->addNeighbour(b->getNeighbour(this), 
                        yhigh,
                        0, end - nx,
                        b->getShift(yhigh)+nx, b->zshift);
      }
    }
  }
  
  // Add boundary between domains
  addNeighbour(d, xhigh, 0, ny-1);
  
  return d;
}

Domain* Domain::splitY(int yind) {
  
  assert( (yind > 0) && (yind < ny) );

  cout << "SplitY " << this << " (" << nx << "," << ny << ") " << yind << "\n";

  Domain *d = create(nx, yind);
  ny -= yind;
  
  cout << "   " << this << " (" << nx << "," << ny << ")\n";
  cout << "   " << d << " (" << d->nx << "," << d->ny << ")\n";

  // Copy the list of boundaries
  list<Bndry*> oldboundary(boundary);
  
  for(list<Bndry*>::iterator it=oldboundary.begin(); it != oldboundary.end(); it++) {
    Bndry *b = *it;
    if(b->onSide(this, yhigh)) {
      // Move to new domain
      b->setNeighbour(ylow, d);
    }else if(b->onSide(this, xlow)) {
      int start = b->getMin(xlow);
      int end   = b->getMax(xlow);
      
      if(start >= ny) {
        // Move to new domain
        b->shiftInds(xlow, -ny);
        b->setNeighbour(xhigh, d);
      }else if(end >= ny) {
        // Crosses cut, so need to split into two
        Domain* neighbour = b->getNeighbour(this);

        b->end -= end - ny + 1;  // Make boundary smaller
        
        // Create a new boundary
        d->addNeighbour(neighbour, 
                        xlow,
                        0, end-ny,
                        b->getShift(xlow)+ny, b->zshift);

      }
    }else if(b->onSide(this, xhigh)) {
      int start = b->getMin(xhigh);
      int end   = b->getMax(xhigh);
      
      if(start >= ny) {
        // Move to new domain
        b->shiftInds(xhigh, -ny);
        b->setNeighbour(xlow, d);
      }else if(end >= ny) {
        // Crosses cut, so need to split into two
        
        b->end -= end - ny + 1; // Make boundary smaller

        // Create a new boundary
        d->addNeighbour(b->getNeighbour(this), 
                        xhigh,
                        0, end - ny,
                        b->getShift(xhigh)+ny, b->zshift);
      }
    }
  }
  
  // Add boundary between domains
  addNeighbour(d, yhigh, 0, nx-1);
  
  return d;
}

std::ostream& operator<<(std::ostream &os, const Domain &d) {
  os << "Domain " << &d << endl;
  os << "Size: " << d.nx << " x " << d.ny << std::endl;
  os << "Connections: ";
  if(d.boundary.empty()) {
    os << "None" << std::endl;
  }else {
    os << std::endl;
    for(int i=0;i<4;i++) {
      os << "\tBoundary "<< i << std::endl;
      for(list<Domain::Bndry*>::const_iterator it=d.boundary.begin(); it != d.boundary.end(); it++) {
        Domain::Bndry *b = *it;
        if(b->onSide(&d, (Domain::BndrySide) i)) {
          os << "\t\tRange: " << b->getMin((Domain::BndrySide)i) << " -> " << b->getMax((Domain::BndrySide)i)
             << " to domain " << b->getNeighbour(&d) << std::endl;
        }
      }
    }
  }
  return os;
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
  
  std::cout << d;
  
  d.addNeighbour(&d, Domain::ylow, 0, 9);
  
  std::cout << d;

  Domain *e = d.splitX(5);
  std::cout << d << *e;

  return 0;
}

#endif
