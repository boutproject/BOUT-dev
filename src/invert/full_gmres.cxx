/*!
 * \file full_gmres.cpp
 *
 * \brief Global inversion using GMRES
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

BoutReal norm_vector(const Field3D &b)
{
  
}

BoutReal dot_product(const Field3D &a, const Field3D &b)
{
  
}

void Update(Field3D &x, int it, BoutReal **h, BoutReal *s, BoutReal *y, Field3D *v)
{
  
}

int full_gmres(const Field3D &b, fgfunc A, Field3D &x, void *extra, int m)
{
  Field3D r, w;
  static Field3D *v = NULL;
  BoutReal normb, beta, resid;
  int it;

  if(v == NULL) {
    v = new Field3D[m];
  }

  normb = norm_vector(b);
  if(normb == 0.0)
    normb = 1.0;

  r = b - A(x,extra);

  beta = norm_vector(r);

  if((resid = beta / normb) <= tol) {
    return(0);
  }

  it = 1;

  while(it < itmax) {
    v[0] = r / beta;
    
    s[0] = beta;

    for(itt=0; (itt < m) && (it <= itmax); itt++, it++) {
      w = A(v[itt], extra);

      
    }
  }
}

