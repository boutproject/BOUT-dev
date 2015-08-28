/**************************************************************************
 * Perpendicular Laplacian inversion. Serial code using FFT
 * and tridiagonal solver.
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

class LaplaceSerialTri;

#ifndef __SERIAL_TRI_H__
#define __SERIAL_TRI_H__

#include <invert_laplace.hxx>
#include <dcomplex.hxx>
#include <options.hxx>

class LaplaceSerialTri : public Laplacian {
public:
  LaplaceSerialTri(Options *opt=NULL);
  ~LaplaceSerialTri();

  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C = val; }
  void setCoefD(const Field2D &val) { D = val; }
  void setCoefEx(const Field2D &val) { bout_error("LaplaceSPT does not have Ex coefficient"); }
  void setCoefEz(const Field2D &val) { bout_error("LaplaceSPT does not have Ez coefficient"); }

  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
private:
  // The coefficents in
  // D*grad_perp^2(x) + (1/C)*grad_perp(C*grad_perp(x)) + A*x = b
  Field2D A, C, D;

  // bk   = The fourier transformed of b, where b is one of the inputs in
  //        LaplaceSerialTri::solve()
  // bk1d = The 1d array of bk
  // xk   = The fourier transformed of x, where x the output of
  //        LaplaceSerialTri::solve()
  // xk1d = The 1d array of xk
  dcomplex **bk, *bk1d;
  dcomplex **xk, *xk1d;

  // Coefficents in the tridiagonal solver matrix
  // Following the notation in "Numerical recipes"
  // avec is the lower diagonal of the matrix
  // bvec is the diagonal of the matrix
  // cvec is the upper diagonal of the matrix
  // NOTE: Do not confuse avec, bvec and cvec with the A, C, and D coefficients
  //       above
  dcomplex *avec, *bvec, *cvec;
};

#endif // __SERIAL_TRI_H__
