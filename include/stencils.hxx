/*!************************************************************************
 * \file stencils.hxx
 * 
 * Sets stencils for differencing
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

#ifndef __STENCILS_H__
#define __STENCILS_H__

#include "bout_types.hxx"

/// Defines a set of values in 1D in the neighbourhood of an index
/// Used for calculating derivatives
struct stencil {
  /// stencil 2 each side of the centre -- in effect means M?G > 2 is not supported
  BoutReal mm = BoutNaN, m = BoutNaN, c = BoutNaN, p = BoutNaN, pp = BoutNaN; 
};


template<DIRECTION direction, STAGGER stagger = STAGGER::None, int nGuard = 1, typename FieldType>
void inline populateStencil(stencil &s, const FieldType& f, const typename FieldType::ind_type i){
  static_assert(nGuard == 1 || nGuard == 2,
		"populateStencil currently only supports one or two guard cells"
		);
  
  switch(stagger) {
  case(STAGGER::None):
    if (nGuard == 2) {
      if (direction == DIRECTION::YOrthogonal) {
        s.mm = f.ynext(-2)[i.template minus<2, direction>()];
      } else {
        s.mm = f[i.template minus<2, direction>()];
      }
    }
    if (direction == DIRECTION::YOrthogonal) {
      s.m = f.ynext(-1)[i.template minus<1, direction>()];
    } else {
      s.m = f[i.template minus<1, direction>()];
    }
    s.c = f[i];
    if (direction == DIRECTION::YOrthogonal) {
      s.p = f.ynext(1)[i.template plus<1, direction>()];
    } else {
      s.p = f[i.template plus<1, direction>()];
    }
    if (nGuard == 2) {
      if (direction == DIRECTION::YOrthogonal) {
        s.pp = f.ynext(2)[i.template plus<2, direction>()];
      } else {
        s.pp = f[i.template plus<2, direction>()];
      }
    }
    break;
  case(STAGGER::C2L):
    if (nGuard == 2) {
      if (direction == DIRECTION::YOrthogonal) {
        s.mm = f.ynext(-2)[i.template minus<2, direction>()];
      } else {
        s.mm = f[i.template minus<2, direction>()];
      }
    }
    if (direction == DIRECTION::YOrthogonal) {
      s.m = f.ynext(-1)[i.template minus<1, direction>()];
    } else {
      s.m = f[i.template minus<1, direction>()];
    }
    s.c = f[i];
    s.p = s.c;
    if (direction == DIRECTION::YOrthogonal) {
      s.pp = f.ynext(1)[i.template plus<1, direction>()];
    } else {
      s.pp = f[i.template plus<1, direction>()];
    }
    break;
  case(STAGGER::L2C):
    if (direction == DIRECTION::YOrthogonal) {
      s.mm = f.ynext(-1)[i.template minus<1, direction>()];
    } else {
      s.mm = f[i.template minus<1, direction>()];
    }
    s.m = f[i];
    s.c = s.m;
    if (direction == DIRECTION::YOrthogonal) {
      s.p = f.ynext(1)[i.template plus<1, direction>()];
    } else {
      s.p = f[i.template plus<1, direction>()];
    }
    if (nGuard == 2) {
      if (direction == DIRECTION::YOrthogonal) {
        s.pp = f.ynext(2)[i.template plus<2, direction>()];
      } else {
        s.pp = f[i.template plus<2, direction>()];
      }
    }
    break;
  }
  return;
}

template<DIRECTION direction, STAGGER stagger = STAGGER::None, int nGuard = 1, typename FieldType >
stencil inline populateStencil(const FieldType& f, const typename FieldType::ind_type i){
  stencil s;
  populateStencil<direction, stagger, nGuard, FieldType>(s, f, i);
  return s;
}
#endif /* __STENCILS_H__ */
