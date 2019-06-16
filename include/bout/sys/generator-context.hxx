#pragma once

#include "bout/region.hxx"
#include "bout/assert.hxx"

class BoundaryRegion;
class Mesh;

#include <map>
#include <string>

class Context {
public:
  /// Set using an index. Can be Ind2D, Ind3D or IndPerp.
  template<IND_TYPE N>
  Context(const SpecificInd<N>& i, CELL_LOC loc, Mesh* msh, BoutReal t)
      : Context(i.x(), i.y(), i.z(), loc, msh, t) {}

  // If Ind2D, z should be 0.0 even if ZLOW (Is this sensible?)
  Context(const Ind2D& i, CELL_LOC loc, Mesh* msh, BoutReal t)
    : Context(i.x(), i.y(), 0, (loc == CELL_ZLOW) ? CELL_CENTRE : loc, msh, t) {}
  
  /// Specify a cell index, together with the cell location, mesh and time
  /// 
  Context(int ix, int iy, int iz, CELL_LOC loc, Mesh* msh, BoutReal t);

  /// If constructed without parameters, contains no values (null).
  /// Requesting x,y,z or t should throw an exception
  ///
  /// NOTE: For backward compatibility, all locations are set to zero.
  /// This should be changed in a future release.
  Context() { parameters["x"] = parameters["y"] = parameters["z"] = parameters["t"] = 0.0; };

  /// The location on the boundary
  Context(const BoundaryRegion* bndry, int iz, CELL_LOC loc, BoutReal t, Mesh* msh);
  Context(const BoundaryRegion* bndry, CELL_LOC loc, BoutReal t, Mesh* msh)
      : Context(bndry, 0, loc, t, msh){};

  Context(const Context& other)
    : localmesh(other.localmesh), parameters(other.parameters) {}

  BoutReal x() const {return get("x");}
  BoutReal y() const {return get("y");}
  BoutReal z() const {return get("z");}
  BoutReal t() const {return get("t");}
  
  /// Set the value of a parameter with given name
  Context& set(const std::string& name, BoutReal value) {
    parameters[name] = value;
    return *this;
  }

  /// Set multiple values, by passing alternating strings and values
  ///
  /// eg. set("x", 1, "y", 2)
  template<typename... Args>
  Context& set(const std::string& name, BoutReal value, Args... args) {
    set(name, value);
    return set(args...);
  }
  
  /// Retrieve a value previously set
  BoutReal get(const std::string& name) const {
    return parameters.at(name);
  }

  /// Get the mesh for this context (position)
  /// If the mesh is null this will throw a BoutException (if CHECK >= 1)
  Mesh* getMesh() const {
    ASSERT1(localmesh)
    return localmesh;
  }

private:
  Mesh *localmesh{nullptr}; ///< The mesh on which the position is defined
  
  /// Contains user-set values which can be set and retrieved
  std::map<std::string, BoutReal> parameters;
};
