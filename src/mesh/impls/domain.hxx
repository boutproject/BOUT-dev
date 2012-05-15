
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

#include <iostream>

using std::string;
using std::list;
using std::find;

#include <uncopyable.hxx>

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

  /// Output info to streams. Mainly useful for debugging
  friend std::ostream& operator<<(std::ostream &os, const Domain &d);
  friend class DomainIterator;
  friend class ConstDomainIterator;
private:
  Domain(); ///< declared, not defined
  
  int nx, ny; ///< Size of the domain
  string name; ///< Name or label for this domain
  
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
    int getZshift(const BndrySide &s) const {
      return (s == side) ? zshift : -zshift;
    }
  };
  
  void addBoundary(Bndry *b);
  
  void removeBoundary(Bndry *b);
  void removeBoundary(std::list<Bndry*>::iterator &it);
  
  list<Bndry*> boundary;

  static BndrySide reverse(const BndrySide &side);
  
  
};

#endif // __DOMAIN_H__
