/**************************************************************************
 * Perpendicular Laplacian inversion
 * 
 * This code uses the Parallel Diagonally Dominant algorithm. This is very efficient
 * (constant number of communications), but achieves this by neglecting "small"
 * corrections. For ELM simulations these seem to be non-negligable, hence:
 *
 * CHECK IF THIS ALGORITHM PRODUCES REASONABLE RESULTS FOR YOUR PROBLEM
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

class LaplacePDD;

#ifndef __LAPLACE_PDD_H__
#define __LAPLACE_PDD_H__

#include <invert_laplace.hxx>
#include <options.hxx>

class LaplacePDD : public Laplacian {
public:
  LaplacePDD(Options *opt = NULL) : Laplacian(opt), Acoef(0.0), Ccoef(1.0), Dcoef(1.0), PDD_COMM_XV(123), PDD_COMM_Y(456) {}
  ~LaplacePDD() {}

  using Laplacian::setCoefA;
  void setCoefA(const Field2D &val) final { Acoef = val; }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) final { Ccoef = val; }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D &val) final { Dcoef = val; }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) final {
    throw BoutException("LaplacePDD does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) final {
    throw BoutException("LaplacePDD does not have Ez coefficient");
  }

  using Laplacian::solve;
  const FieldPerp solve(const FieldPerp &b);
  const Field3D solve(const Field3D &b);
private:
  Field2D Acoef, Ccoef, Dcoef;
  
  const int PDD_COMM_XV; // First message tag
  const int PDD_COMM_Y;  // Second tag
  
  /// Data structure for PDD algorithm
  typedef struct {
    dcomplex **bk;  ///< b vector in Fourier space

    dcomplex **avec, **bvec, **cvec; ///< Diagonal bands of matrix
  
    int jy; ///< Y index
  
    dcomplex **xk;
    dcomplex **v, **w;

    BoutReal *snd; // send buffer
    BoutReal *rcv; // receive buffer
  
    comm_handle recv_handle;

    dcomplex *y2i;
  }PDD_data;
  
  void start(const FieldPerp &b, PDD_data &data);
  void next(PDD_data &data);
  void finish(PDD_data &data, FieldPerp &x);
};

#endif // __LAPLACE_PDD_H__
