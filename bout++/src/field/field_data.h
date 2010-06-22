/*!
 * \file field_data.h
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

#ifndef __FIELD_DATA_H__
#define __FIELD_DATA_H__

#include "bout_types.h"

/// Interface used to access data in field classes
/*!
  Used by communicator, solver and (soon) datafile classes
  to access internal data in a general way
*/
class FieldData {
 public:
  virtual ~FieldData() { }

  // Defines interface which must be implemented
  virtual bool isReal() const = 0; ///< Returns true if field consists of real values
  virtual bool is3D() const = 0;   ///< True if variable is 3D
  
  virtual int byteSize() const = 0; ///< Number of bytes for a single point
  virtual int realSize() const { return 0; } ///< Number of reals (not implemented if not real)

  virtual int getData(int x, int y, int z, void *vptr) const = 0; ///< Return number of bytes
  virtual int getData(int x, int y, int z, real *rptr) const = 0; ///< Return number of reals
  
  virtual int setData(int x, int y, int z, void *vptr) = 0;
  virtual int setData(int x, int y, int z, real *rptr) = 0;

  // This code for inputting/outputting to file (all optional)
  virtual bool  ioSupport() { return false; }  ///< Return true if these functions implemented
  virtual const string getSuffix(int component) const { return string(""); }
  virtual void* getMark() const {return NULL;} ///< Store current settings (e.g. co/contra-variant)
  virtual void  setMark(void *setting) {}      ///< Return to the stored settings
  virtual real* getData(int component) { return NULL; }
  virtual void  zeroComponent(int component) { } ///< Set a component to zero
  
  /// Added 20/8/2008 for twist-shifting in communication routine
  virtual void shiftZ(int jx, int jy, double zangle) { }

#ifdef CHECK
  virtual void doneComms() { }; // Notifies that communications done
#endif

};

#endif

