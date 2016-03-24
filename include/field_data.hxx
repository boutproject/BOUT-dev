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
#include <string>
using std::string;

// Including the next line leads to compiler errors
//#include "boundary_op.hxx"
class BoundaryOp;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <boundary_region.hxx>

class FieldGenerator; // Forward declaration

/// Interface used to access data in field classes
/*!
  Used by communicator, solver and (soon) datafile classes
  to access internal data in a general way
*/
class FieldData {
 public:
  FieldData();
  virtual ~FieldData();

  // Defines interface which must be implemented
  virtual bool isReal() const = 0; ///< Returns true if field consists of BoutReal values
  virtual bool is3D() const = 0;   ///< True if variable is 3D
  
  virtual int byteSize() const = 0; ///< Number of bytes for a single point
  virtual int BoutRealSize() const { return 0; } ///< Number of BoutReals (not implemented if not BoutReal)

  virtual int getData(int x, int y, int z, void *vptr) const = 0; ///< Return number of bytes
  virtual int getData(int x, int y, int z, BoutReal *rptr) const = 0; ///< Return number of BoutReals
  
  virtual int setData(int x, int y, int z, void *vptr) = 0;
  virtual int setData(int x, int y, int z, BoutReal *rptr) = 0;
  
  /// Added 20/8/2008 for twist-shifting in communication routine
  virtual void shiftZ(int jx, int jy, double zangle) { }

#ifdef CHECK
  virtual void doneComms() { }; // Notifies that communications done
#endif
  
  // Boundary conditions
  void setBoundary(const string &name); ///< Set the boundary conditions
  void setBoundary(const string &region, BoundaryOp *op); ///< Manually set

  void copyBoundary(const FieldData &f); ///< Copy the boundary conditions from another field

  virtual void applyBoundary() {}
  virtual void applyTDerivBoundary() {};
//JMAD
  void addBndryFunction(FuncPtr userfunc, BndryLoc location);
  void addBndryGenerator(FieldGenerator* gen, BndryLoc location);
  
  FieldGenerator* getBndryGenerator(BndryLoc location);

 protected:
  vector<BoundaryOp*> bndry_op; // Boundary conditions
  bool boundaryIsCopy; // True if bndry_op is a copy
  
  bool boundaryIsSet; // Set to true when setBoundary called
 
 std::map <BndryLoc,FieldGenerator*> bndry_generator;
};

#endif

