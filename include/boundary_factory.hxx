
class BoundaryFactory;

#ifndef __BNDRY_FACTORY_H__
#define __BNDRY_FACTORY_H__

#include "boundary_op.hxx"
#include "boundary_region.hxx"
#include "parallel_boundary_op.hxx"
#include "parallel_boundary_region.hxx"

#include <string>
#include <map>
using std::string;
using std::map;

/// Create BoundaryOp objects on demand
class BoundaryFactory {
 public:
  ~BoundaryFactory();
  /// Return a pointer to the only instance
  static BoundaryFactory* getInstance();

  static void cleanup(); ///< Frees all memory

  /// Create a boundary operation object
  BoundaryOpBase* create(const string &name, BoundaryRegionBase *region);
  BoundaryOpBase* create(const char* name, BoundaryRegionBase *region);

  /// Create a boundary object using the options file
  BoundaryOpBase* createFromOptions(const string &varname, BoundaryRegionBase *region);
  BoundaryOpBase* createFromOptions(const char* varname, BoundaryRegionBase *region);

  // functions to add available boundary conditions and modifiers
  // Supply an object, the name, and (optionally) which boundaries it can be applied to
  void add(BoundaryOp* bop, const string &name);
  void add(BoundaryOp* bop, const char *name);
  void addMod(BoundaryModifier* bmod, const string &name);
  void addMod(BoundaryModifier* bmod, const char *name);
  // Parallel boundaries
  void add(BoundaryOpPar* bop, const string &name);
  void add(BoundaryOpPar* bop, const char *name);

 private:
  BoundaryFactory(); // Prevent instantiation of this class
  static BoundaryFactory* instance; ///< The only instance of this class (Singleton)

  // Database of available boundary conditions and modifiers
  map<string, BoundaryOp*> opmap;
  map<string, BoundaryModifier*> modmap;
  // Parallel boundary conditionscccc
  map<string, BoundaryOpPar*> par_opmap;
  // Modifiers to be implemented...
  // map<string, BoundaryModifier*> par_modmap;

  // Functions to look up operations and modifiers
  BoundaryOp* findBoundaryOp(const string &s);
  BoundaryModifier* findBoundaryMod(const string &s);
  // Parallel boundary conditions
  BoundaryOpPar* findBoundaryOpPar(const string &s);
  // To be implemented...
  // BoundaryModifier* findBoundaryMod(const string &s);

};

#endif // __BNDRY_FACTORY_H__

