
class Domain;

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include <list>

#include <uncopyable.hxx>

/// Logically rectangular domain
class Domain : private Uncopyable {
public:
  Domain(int NX, int NY);
  ~Domain();
  
  /// Boundary locations
  enum BndrySide { xlow=0, xhigh=1, ylow=2, yhigh=3 };
  
  /// Add a neighbour to the domain
  void addNeighbour(Domain *d, BndrySide side, int first, int last, int shift=0);

  Domain* splitX(int xind);
  Domain* splitY(int yind);
  
private:
  Domain(); ///< declared, not defined
  
  int nx, ny; ///< Size of the domain
  
  struct Bndry {
    BndrySide side; ///< Side of the domain
    int start, end;  ///< Index range
    Domain *neighbour;
  };
  
  void removeNeighbour(Domain *d);
  void removeNeighbour(Domain *d, BndrySide side);
  void removeNeighbour(std::list<Bndry>::iterator &it);

  std::list<Bndry> boundary;

  static BndrySide reverse(const BndrySide &side);
};

#endif // __DOMAIN_H__
