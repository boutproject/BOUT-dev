
class BoundaryFactory;

#ifndef __BNDRY_FACTORY_H__
#define __BNDRY_FACTORY_H__

#include "boundary_op.hxx"
#include "boundary_region.hxx"
#include "parallel_boundary_op.hxx"
#include "parallel_boundary_region.hxx"

#include <string>
#include <map>

/// Create BoundaryOp objects on demand
/*!
 * This implements a simple string parser, used to match boundary condition
 * names like "dirichlet" with a BoundaryOp object.
 * 
 * Modifiers: Simple modifications of boundary conditions can be performed,
 * for example transforming the coordinate system
 * 
 * This is a singleton, so only one instance can exist. This is 
 * enforced by making the constructor private, and having a getInstance() method
 * to return a pointer to the only instance.
 *
 * Example
 * -------
 *
 * Boundaries are defined as classes which inherit from BoundaryOp
 * These define a clone() function which creates a new BoundaryOp,
 * given a list of arguments. See boundary_standard.hxx for examples.
 * 
 *     class MyBoundary : public BoundaryOp {
 *      public:
 *       BoundaryOp* clone(BoundaryRegion *region, const list<string> &args) {
 *         // Decide what to do with arguments
 *         return new MyBoundary();
 *       }
 *       void apply(Field2D &f);
 *       void apply(Field3D &f);
 *     };
 * 
 * The singleton instance of BoundaryFactory from getInstance():
 * 
 *     BoundaryFactory* bf = BoundaryFactory::getInstance();
 *
 * New boundary types can be added to the BoundaryFactory
 *
 *     bf->add(new MyBoundary, "myboundary");
 * 
 *
 * Subsequent calls to create() or createFromOptions() can make use
 * of the boundary type "myboundary".
 *
 * BoundaryOpBase *bndry = bf->create("myboundary()", new BoundaryRegionXOut("xout", 0, 10, localmesh));
 * 
 * where the region is defined in boundary_region.hxx
 * 
 */
class BoundaryFactory {
 public:
  ~BoundaryFactory();
  /// Return a pointer to the only instance
  static BoundaryFactory* getInstance();

  static void cleanup(); ///< Frees all memory

  /// Create a boundary operation object
  BoundaryOpBase* create(const std::string &name, BoundaryRegionBase *region);
  BoundaryOpBase* create(const char* name, BoundaryRegionBase *region);

  /// Create a boundary object using the options file
  BoundaryOpBase* createFromOptions(const std::string &varname, BoundaryRegionBase *region);
  BoundaryOpBase* createFromOptions(const char* varname, BoundaryRegionBase *region);

  /*!
   * Add available boundary conditions and modifiers
   * Supply an object, and the name to be used
   */
  void add(BoundaryOp* bop, const std::string &name);
  
  /*!
   * Add a boundary condition.
   *
   * Note: This method should be removed, as the string method is sufficient
   */
  void add(BoundaryOp* bop, const char *name);
  
  /*!
   * Add a boundary condition modifier
   */
  void addMod(BoundaryModifier* bmod, const std::string &name);
  
  /*!
   * Note: This method should be removed, as the string method is sufficient
   */ 
  void addMod(BoundaryModifier* bmod, const char *name);
  
  // Parallel boundaries
  void add(BoundaryOpPar* bop, const std::string &name);
  void add(BoundaryOpPar* bop, const char *name);

 private:
  /*!
   * Private constructor, preventing instantiation of this class
   */ 
  BoundaryFactory(); 
  static BoundaryFactory* instance; ///< The only instance of this class (Singleton)

  // Database of available boundary conditions and modifiers
  std::map<std::string, BoundaryOp*> opmap;
  std::map<std::string, BoundaryModifier*> modmap;
  // Parallel boundary conditions
  std::map<std::string, BoundaryOpPar*> par_opmap;
  // Modifiers to be implemented...
  // map<string, BoundaryModifier*> par_modmap;

  // Functions to look up operations and modifiers
  BoundaryOp* findBoundaryOp(const std::string &s);
  BoundaryModifier* findBoundaryMod(const std::string &s);
  // Parallel boundary conditions
  BoundaryOpPar* findBoundaryOpPar(const std::string &s);
  // To be implemented...
  // BoundaryModifier* findBoundaryMod(const string &s);

};

#endif // __BNDRY_FACTORY_H__

