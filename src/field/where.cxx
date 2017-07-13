/**************************************************************************
 * A set of functions which choose between two values
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

#include <where.hxx>

const Field3D where(const Field2D &test, const Field3D &gt0, const Field3D &le0) {
  Field3D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0[i];
    }else {
      result[i] = le0[i];
    }
  }
  return result;
}

const Field3D where(const Field2D &test, const Field3D &gt0, BoutReal le0) {
  Field3D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0[i];
    }else {
      result[i] = le0;
    }
  }
  return result;
}

const Field3D where(const Field2D &test, BoutReal gt0, const Field3D &le0) {
  Field3D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0;
    }else {
      result[i] = le0[i];
    }
  }
  
  return result;
}

const Field3D where(const Field2D &test, const Field3D &gt0, const Field2D &le0) {
  Field3D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0[i];
    }else {
      result[i] = le0[i];
    }
  }
  
  return result;
}

const Field3D where(const Field2D &test, const Field2D &gt0, const Field3D &le0) {
  Field3D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0[i];
    }else {
      result[i] = le0[i];
    }
  }
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////////
// Versions returning Field2D

const Field2D where(const Field2D &test, const Field2D &gt0, const Field2D &le0) {
  Field2D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0[i];
    }else {
      result[i] = le0[i];
    }
  }
  
  return result;
}

const Field2D where(const Field2D &test, const Field2D &gt0, BoutReal le0) {
  Field2D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0[i];
    }else {
      result[i] = le0;
    }
  }
  
  return result;
}

const Field2D where(const Field2D &test, BoutReal gt0, const Field2D &le0) {
  Field2D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0;
    }else {
      result[i] = le0[i];
    }
  }
  
  return result;
}

const Field2D where(const Field2D &test, BoutReal gt0, BoutReal le0) {
  Field2D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0;
    }else {
      result[i] = le0;
    }
  }
  
  return result;
}

const Field3D where(const Field3D &test, BoutReal gt0, const Field3D &le0) {
  Field3D result(test.getMesh());
  result.allocate();
  
  for(auto i : result) {
    if(test[i] > 0.0) {
      result[i] = gt0;
    }else {
      result[i] = le0[i];
    }
  }
  return result;
}
