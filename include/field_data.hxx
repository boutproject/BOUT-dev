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
#ifndef __FIELD_DATA_H__
#define __FIELD_DATA_H__

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

#include "boundary_region.hxx"
#include "parallel_boundary_region.hxx"

#include "bout/sys/expressionparser.hxx"

/// Interface used to access data in field classes
/*!
  Used by communicator, solver and (soon) datafile classes
  to access internal data in a general way
*/
class FieldData {
public:
  FieldData() = default;
  FieldData(const FieldData& other);
  FieldData(FieldData&& other) = default;
  FieldData& operator=(const FieldData& other);
  FieldData& operator=(FieldData&& other) = default;
  virtual ~FieldData();

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

protected:
  std::vector<BoundaryOp *> bndry_op; ///< Boundary conditions
  bool boundaryIsCopy{false};         ///< True if bndry_op is a copy
  bool boundaryIsSet{false};          ///< Set to true when setBoundary called

  // Parallel boundaries
  std::vector<BoundaryOpPar *> bndry_op_par; ///< Boundary conditions

  std::map <BndryLoc,FieldGeneratorPtr> bndry_generator;
};

#endif

