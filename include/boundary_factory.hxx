
class BoundaryFactory;

#ifndef __BNDRY_FACTORY_H__
#define __BNDRY_FACTORY_H__

#include "boundary_op.hxx"
#include "boundary_region.hxx"

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
  BoundaryOp* create(const string &name, BoundaryRegion *region);
  BoundaryOp* create(const char* name, BoundaryRegion *region);

  /// Create a boundary object using the options file
  BoundaryOp* createFromOptions(const string &varname, BoundaryRegion *region);
  BoundaryOp* createFromOptions(const char* varname, BoundaryRegion *region);

  // functions to add available boundary conditions and modifiers
  // Supply an object, the name, and (optionally) which boundaries it can be applied to
  void add(BoundaryOp* bop, const string &name);
  void add(BoundaryOp* bop, const char *name);
  void addMod(BoundaryModifier* bmod, const string &name);
  void addMod(BoundaryModifier* bmod, const char *name);

 private:
  BoundaryFactory(); // Prevent instantiation of this class
  static BoundaryFactory* instance; ///< The only instance of this class (Singleton)

  // Database of available boundary conditions and modifiers
  map<string, BoundaryOp*> opmap;
  map<string, BoundaryModifier*> modmap;

  // Functions to look up operations and modifiers
  BoundaryOp* findBoundaryOp(const string &s);
  BoundaryModifier* findBoundaryMod(const string &s);
};

#endif // __BNDRY_FACTORY_H__

