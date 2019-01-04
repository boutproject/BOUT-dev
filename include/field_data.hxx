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

#include "boundary_region.hxx" // Just for BndryLoc enum
#include "bout_types.hxx"
#include "unused.hxx"
#include "bout/sys/expressionparser.hxx"

#include <map>
#include <memory>
#include <string>
#include <vector>

class BoundaryOp;
class BoundaryOpPar;
class FieldVisitor;

/// Interface used to access data in field classes
/*!
  Used by communicator, solver and (soon) datafile classes
  to access internal data in a general way
*/
class FieldData {
public:
  FieldData() = default;
  virtual ~FieldData();

  // Visitor pattern support
  virtual void accept(FieldVisitor &v) = 0;
  
  // Defines interface which must be implemented
  virtual bool isReal() const = 0; ///< Returns true if field consists of BoutReal values
  virtual bool is3D() const = 0;   ///< True if variable is 3D
  
  virtual int byteSize() const = 0; ///< Number of bytes for a single point
  virtual int BoutRealSize() const { return 0; } ///< Number of BoutReals (not implemented if not BoutReal)

  virtual void doneComms() { }; // Notifies that communications done
  
  // Boundary conditions
  void setBoundary(const std::string &name); ///< Set the boundary conditions
  void setBoundary(const std::string &region, BoundaryOp *op); ///< Manually set

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

