
#include "domain.hxx"

#include <boutexception.hxx>

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
  
  std::cout << "Adding boundary: " << b << endl;
  
  // Add to both domains
  addBoundary(b);
  if((d != NULL) && (d != this))
    d->addBoundary(b);
}

void Domain::addBoundary(Bndry *b) {
  // Error checking
  
  // Add to list of boundaries
  boundary.push_back(b);
}

void Domain::removeBoundary(Bndry *b) {
  boundary.remove(b);
}

void Domain::removeBoundary(list<Bndry*>::iterator &it) {
  it = boundary.erase(it);
}


Domain* Domain::splitX(int xind) {
  Domain *d = new Domain(xind, ny); // New domain on the left
  
  nx -= xind;
  
  for(list<Bndry*>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    Bndry *b = *it;
    
    if(b->onSide(this, xlow)) {
      // Move to new domain
      b->setNeighbour(xhigh, d);
      it = boundary.begin(); // Start at the beginning again
    }else if(b->onSide(this, ylow)) {
      int start = b->getMin(ylow);
      int end   = b->getMax(ylow);
      
      if(end < xind) {
        // Move to new domain
        b->setNeighbour(yhigh, d);
        it = boundary.begin();
      }else if((start <= xind) &&
               (end >= xind)) {
        // Crosses cut, so need to split into two
        
        b->setNeighbour(yhigh, d); // Move to new domain
        b->end -= end - xind + 1;      // Make boundary smaller
        
        // Create a new boundary
        addNeighbour(b->getNeighbour(this), 
                     ylow,
                     0, end - xind,
                     b->shift+xind, b->zshift);
        it = boundary.begin();
      }
    }else if(b->onSide(this, yhigh)) {
      int start = b->getMin(yhigh);
      int end   = b->getMax(yhigh);
      
      if(end < xind) {
        // Move to new domain
        b->setNeighbour(ylow, d);
        it = boundary.begin();
      }else if((start <= xind) &&
               (end >= xind)) {
        // Crosses cut, so need to split into two
        
        b->setNeighbour(ylow, d); // Move to new domain
        b->end -= end - xind + 1;      // Make boundary smaller
        
        // Create a new boundary
        addNeighbour(b->getNeighbour(this), 
                     yhigh,
                     0, end - xind,
                     b->shift+xind, b->zshift);
        it = boundary.begin();
      }
    }
  }
  
  // Add boundary between domains
  addNeighbour(d, xlow, 0, ny-1);
  
  return d;
}

Domain* Domain::splitY(int yind) {
  Domain *d = new Domain(nx, yind);
  
  ny -= yind;

  for(list<Bndry*>::iterator it=boundary.begin(); it != boundary.end(); it++) {
    Bndry *b = *it;
    cout << b << endl;
    if(b->onSide(this, ylow)) {
      // Move to new domain
      b->setNeighbour(yhigh, d);
      it = boundary.begin();
    }else if(b->onSide(this, xlow)) {
      int start = b->getMin(xlow);
      int end   = b->getMax(xlow);
      
      if(end < yind) {
        // Move to new domain
        b->setNeighbour(xhigh, d);
        it = boundary.begin();
      }else if((start <= yind) &&
               (end >= yind)) {
        // Crosses cut, so need to split into two
        
        b->setNeighbour(xhigh, d); // Move to new domain
        b->end -= end - yind + 1;      // Make boundary smaller
        
        // Create a new boundary
        addNeighbour(b->getNeighbour(this), 
                     xlow,
                     0, end-yind,
                     b->shift+yind, b->zshift);
        it = boundary.begin();
      }
    }else if(b->onSide(this, xhigh)) {
      int start = b->getMin(xhigh);
      int end   = b->getMax(xhigh);
      
      if(end < yind) {
        // Move to new domain
        b->setNeighbour(xlow, d);
        it = boundary.begin();
      }else if((start <= yind) &&
               (end >= yind)) {
        // Crosses cut, so need to split into two
        
        b->setNeighbour(xlow, d); // Move to new domain
        b->end -= end - yind + 1;      // Make boundary smaller
        
        // Create a new boundary
        addNeighbour(b->getNeighbour(this), 
                     xhigh,
                     0, end - yind,
                     b->shift+yind, b->zshift);
        it = boundary.begin();
      }
    }
  }
  
  // Add boundary between domains
  addNeighbour(d, ylow, 0, ny-1);
  

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
