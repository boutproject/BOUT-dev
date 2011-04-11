/**************************************************************************
 * Maintains a list of values which are printed every time-step
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
 **************************************************************************/

#ifndef __DIAGNOS_H__
#define __DIAGNOS_H__

#include <vector>
using std::vector;

#include "field_data.hxx"

/// Functions which can be applied to the data
enum DIAGNOS_FUNC {DIAG_INDX, DIAG_MAX, DIAG_MAXABS, DIAG_MIN, DIAG_MEAN, DIAG_RMS, DIAG_MAX_ZRMS};

class Diagnos {
 public:
  Diagnos();
  ~Diagnos();
  void add(FieldData &f, DIAGNOS_FUNC func, int x=-1, int y=-1, int z=-1, int component=-1, const char *label = NULL);
  
  const vector< BoutReal > run();

 private:

  static bool init;
  static bool global_vals; ///< If true, prints global values

  typedef struct {
    FieldData *var;

    char *label;

    DIAGNOS_FUNC func;
    int x, y, z;
    int component;
  }diag_item;

  vector< diag_item > item;

  BoutReal run(const diag_item &i);
};


#endif // __DIAGNOS_H__
