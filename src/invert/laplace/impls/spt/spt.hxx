/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           PARALLEL CODE - SIMPLE ALGORITHM
 * 
 * I'm just calling this Simple Parallel Tridag. Naive parallelisation of
 * the serial code. For use as a reference case.
 * 
 * Overlap calculation / communication of poloidal slices to achieve some
 * parallelism.
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

class LaplaceSPT;

#ifndef __SPT_H__
#define __SPT_H__

#include <invert_laplace.hxx>
#include <dcomplex.hxx>
#include <options.hxx>

/// Simple parallelisation of the Thomas tridiagonal solver algorithm (serial code)
/*!
 * This is a reference code which performs the same operations as the serial code.
 * To invert a single XZ slice (FieldPerp object), data must pass from the innermost
 * processor (mesh->PE_XIND = 0) to the outermost (mesh->PE_XIND = mesh->NXPE-1) and back again.
 *
 * Some parallelism is achieved by running several inversions simultaneously, so while
 * processor #1 is inverting Y=0, processor #0 is starting on Y=1. This works ok as long
 * as the number of slices to be inverted is greater than the number of X processors (MYSUB > mesh->NXPE).
 * If MYSUB < mesh->NXPE then not all processors can be busy at once, and so efficiency will fall sharply.
 *
 * @param[in]    b      RHS values (Ax = b)
 * @param[in]    flags  Inversion settings (see boundary.h for values)
 * @param[in]    a      This is a 2D matrix which allows solution of A = Delp2 + a
 * @param[out]   data   Structure containing data needed for second half of inversion
 * @param[in]    ccoef  Optional coefficient for first-order derivative
 * @param[in]    d      Optional factor to multiply the Delp2 operator
 */
class LaplaceSPT : public Laplacian {
public:
  LaplaceSPT(Options *opt = NULL) : Laplacian(opt), A(0.0), C(1.0), D(1.0), SPT_DATA(1123) {}
  ~LaplaceSPT() {}
  
  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C = val; }
  void setCoefD(const Field2D &val) { D = val; }
  void setCoefEx(const Field2D &val) { bout_error("LaplaceSPT does not have Ex coefficient"); }
  void setCoefEz(const Field2D &val) { bout_error("LaplaceSPT does not have Ez coefficient"); }
  
  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
  
  const Field3D solve(const Field3D &b);
  const Field3D solve(const Field3D &b, const Field3D &x0);
private:
  Field2D A, C, D;
  
  /// Data structure for SPT algorithm
  typedef struct {
    int jy; ///< Y index
    
    dcomplex **bk;  ///< b vector in Fourier space
    dcomplex **xk;

    dcomplex **gam;
  
    dcomplex **avec, **bvec, **cvec; ///< Diagonal bands of matrix

    int proc; // Which processor has this reached?
    int dir;  // Which direction is it going?
  
    comm_handle recv_handle; // Handle for receives
  
    int comm_tag; // Tag for communication
  
    BoutReal *buffer;
  }SPT_data;
  
  void tridagForward(dcomplex *a, dcomplex *b, dcomplex *c,
                      dcomplex *r, dcomplex *u, int n,
                      dcomplex *gam,
                      dcomplex &bet, dcomplex &um, bool start=false);
  void tridagBack(dcomplex *u, int n,
                   dcomplex *gam, dcomplex &gp, dcomplex &up);
  
  const int SPT_DATA; ///< 'magic' number for SPT MPI messages
  
  int start(const FieldPerp &b, SPT_data &data);
  
  int next(SPT_data &data);
  
  void finish(SPT_data &data, FieldPerp &x);

};

#endif // __SPT_H__
