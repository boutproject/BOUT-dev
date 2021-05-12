/*!
 * \file field_data.hxx
 * \brief Class inherited by any field wanting to use Communicator or Solver objects
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

class FieldData;

#pragma once
#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include "bout_types.hxx"
#include "unused.hxx"

#include <map>
#include <memory>
#include <string>
#include <vector>

// Including the next line leads to compiler errors
//#include "boundary_op.hxx"
class BoundaryOp;
class BoundaryOpPar;
class Coordinates;
class Mesh;

#include "boundary_region.hxx"
#include "parallel_boundary_region.hxx"

#include "bout/sys/expressionparser.hxx"

/// Base class for both scalar and vector fields, holds common
/// information about the grid, coordinates, and boundaries
class FieldData {
public:
  FieldData() = default;
  FieldData(const FieldData& other);
  FieldData(FieldData&& other) = default;
  FieldData& operator=(const FieldData& other);
  FieldData& operator=(FieldData&& other) = default;
  virtual ~FieldData();

  FieldData(Mesh* localmesh, CELL_LOC location_in = CELL_LOC::centre);

  /// Set variable location for staggered grids to @param new_location
  ///
  /// Throws BoutException if new_location is not `CELL_CENTRE` and
  /// staggered grids are turned off and checks are on. If checks are
  /// off, silently sets location to ``CELL_CENTRE`` instead.
  virtual FieldData& setLocation(CELL_LOC new_location);
  /// Get variable location
  virtual CELL_LOC getLocation() const;

  // Defines interface which must be implemented
  /// True if variable is 3D
  virtual bool is3D() const = 0;
  /// Number of BoutReals in one element
  virtual int elementSize() const { return 1; }

  virtual void doneComms() { }; // Notifies that communications done
  
  // Boundary conditions
  void setBoundary(const std::string &name); ///< Set the boundary conditions

  void copyBoundary(const FieldData &f); ///< Copy the boundary conditions from another field

  virtual void applyBoundary(bool UNUSED(init)=false) {}
  virtual void applyTDerivBoundary() {};
//JMAD
  void addBndryFunction(FuncPtr userfunc, BndryLoc location);
  void addBndryGenerator(FieldGeneratorPtr gen, BndryLoc location);
  
  FieldGeneratorPtr getBndryGenerator(BndryLoc location);

  /// Returns a pointer to the `Mesh` object used by this field
  Mesh* getMesh() const;

  /// Returns a pointer to the `Coordinates` object at this field's
  /// location from the mesh this field is on.
  Coordinates* getCoordinates() const;

  /// Returns a pointer to the `Coordinates` object at the requested
  /// location from the mesh this field is on. If \p loc is `CELL_DEFAULT`
  /// then return coordinates at field location
  Coordinates* getCoordinates(CELL_LOC loc) const;

  /// Has `setBoundary` been called?
  bool isBoundarySet() const { return boundaryIsSet; }

  /// Get boundary conditions
  std::vector<BoundaryOp*> getBoundaryOps() const { return bndry_op; }
  /// Get parallel boundary conditions
  std::vector<BoundaryOpPar*> getBoundaryOpPars() const { return bndry_op_par; }

  friend void swap(FieldData& first, FieldData& second) noexcept;

protected:
  /// Grid information, etc. Owned by the simulation or global object
  Mesh* fieldmesh{nullptr};

private:
  std::vector<BoundaryOp *> bndry_op; ///< Boundary conditions
  bool boundaryIsCopy{false};         ///< True if bndry_op is a copy
  bool boundaryIsSet{false};          ///< Set to true when setBoundary called

  // Parallel boundaries
  std::vector<BoundaryOpPar *> bndry_op_par; ///< Boundary conditions

  std::map <BndryLoc,FieldGeneratorPtr> bndry_generator;

  /// `Coordinates` used by this field, owned by `fieldmesh`
  mutable std::weak_ptr<Coordinates> fieldCoordinates{};
  /// Location of the variable in the cell
  CELL_LOC location{CELL_CENTRE};
};

#endif
