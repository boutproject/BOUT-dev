
class Domain;
class DomainIterator;
class ConstDomainIterator;

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include <list>
#include <ostream>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>

using std::string;
using std::list;
using std::find;

#include <bout/sys/uncopyable.hxx>

class DomainIterator : public std::iterator<std::forward_iterator_tag, Domain> {
public:
  DomainIterator();
  DomainIterator(Domain &d);
  
  Domain& operator*() {return *dom;}
  Domain* operator->() {return dom;}

  bool operator==(const DomainIterator &x) const {
    return dom == x.dom;
  }
  bool operator!=(const DomainIterator &x) const {
    return dom != x.dom;
  }
  DomainIterator& operator++();
  DomainIterator operator++(int) {DomainIterator tmp(*this); operator++(); return tmp;}
private:
  Domain* dom; // Current domain
  list<Domain*> visited_list; // Previously visited (includes dom)
  list<Domain*> stack; // Stack of visited domains (includes dom)
  
  bool visited(Domain *d) {
    list<Domain*>::iterator it = find(visited_list.begin(), visited_list.end(), d);
    return it != visited_list.end();
  }
};

class ConstDomainIterator {
public:
  ConstDomainIterator() : dom(0) {}
  ConstDomainIterator(const Domain &d) : dom(&d) {}
  
  const Domain& operator*() const {return *dom;}
  Domain const* operator->() const {return dom;}
  
  bool operator==(const ConstDomainIterator &x) const {
    return dom == x.dom;
  }
  bool operator!=(const ConstDomainIterator &x) const {
    return dom != x.dom;
  }
  
  ConstDomainIterator& operator++();
private:
  const Domain* dom;
};

/// Logically rectangular domain
class Domain : private Uncopyable {
public:
  Domain(int NX, int NY);
  Domain(int NX, int NY, const string &Name);
  ~Domain();

  // Polymorphic constructors
  virtual Domain* create(int NX, int NY) { return new Domain(NX, NY); }
  
  // Iterators over domains
  typedef DomainIterator iterator;
  typedef ConstDomainIterator const_iterator;

  iterator  begin() {return iterator(*this);}
  const_iterator begin() const {return const_iterator(*this);}
  
  iterator end() {return iterator();}
  const_iterator end() const {return const_iterator();}
  
  /// Boundary locations
  enum BndrySide { xlow=0, xhigh=1, ylow=2, yhigh=3 };
  
  /// Add a neighbour to the domain
  void addNeighbour(Domain *d, BndrySide side, int first, int last, int shift=0, int zshift=0);

  /// Split domains in two
  Domain* splitX(int xind);
  Domain* splitY(int yind);
  
  int xSize() const {return nx;}
  int ySize() const {return ny;}
  int area() const {return nx*ny;}

  // Allow user to attach extra 
  void setUserData(void* data) {user_data = data;}
  void* getUserData() const {return user_data;}

  // Keep track of a global origin
  void setOrigin(int x, int y) {x0 = x; y0 = y;}
  int xOrigin() const {return x0;}
  int yOrigin() const {return y0;}
  
  /// Output info to streams. Mainly useful for debugging
  friend std::ostream& operator<<(std::ostream &os, const Domain &d);
  friend class DomainIterator;
  friend class ConstDomainIterator;
  
  /// Boundary structure, shared between domains
  struct Bndry {
  public:
    Bndry(BndrySide s, Domain *d, Domain *other, int f, int l, int sh=0, int zsh=0)
      : side(s), from(d), to(other), start(f), end(l), shift(sh), zshift(zsh) {}
    
    BndrySide side;  ///< Side of the domain 'from'
    Domain *from, *to; ///< Domains either side of the boundary
    int start, end;  ///< Index range in 'from' domain
    int shift;       ///< Shift going from this domain to neighbour
    int zshift;      ///< Shift in Z across boundary

    // Check if boundary is on a given side
    // Done this way as boundary could be on two sides if periodic
    bool onSide(const Domain *me, BndrySide s) const {
      if( (me == from) && (s == side) )
        return true;
      if( (me == to) && (s == reverse(side)) )
        return true;
      return false;
    }
    int getMin(BndrySide s) {
      return (s == side) ? start : start+shift;
    }
    int getMax(BndrySide s) {
      return (s == side) ? end : end + shift;
    }

    const Domain* getNeighbour(const Domain* me) const {
      return (me == from) ? to : from;
    }
    Domain* getNeighbour(Domain* me) const {
      return (me == from) ? to : from;
    }
    void setNeighbour(BndrySide s, Domain *other) {
      if(s == side) {
        if((to != 0) && (from != to))
          to->removeBoundary(this);
        to = other;
        if(to != 0)
          to->addBoundary(this);
      }else {
        if((from != 0) && (from != to))
          from->removeBoundary(this);
        from = other;
        if(from != 0)
          from->addBoundary(this);
      }
    }

    // Shift the indices in only one of the domains
    void shiftInds(BndrySide s, int x) {
      if(s == side) {
        start += x;
        end += x;
        shift -= x;
      }else {
        shift += x;
      }
    }
    
    int getShift(BndrySide s) {
      return (s == side) ? shift : -shift;
    }
    
    int getZshift(const BndrySide &s) const {
      return (s == side) ? zshift : -zshift;
    }
  };
  
  /// Iterator over boundaries
  typedef list<Bndry*>::const_iterator bndry_iterator;
  bndry_iterator bndry_begin() const { return boundary.begin(); }
  bndry_iterator bndry_end() const { return boundary.end(); }
private:
  Domain(); ///< declared, not defined
  
  int x0, y0; ///< Origin (for user)
  int nx, ny; ///< Size of the domain
  string name; ///< Name or label for this domain
  
  void addBoundary(Bndry *b);
  void removeBoundary(Bndry *b);
  void removeBoundary(std::list<Bndry*>::iterator &it);
  
  list<Bndry*> boundary;

  static BndrySide reverse(const BndrySide &side);
  
  void* user_data;
};

#endif // __DOMAIN_H__
