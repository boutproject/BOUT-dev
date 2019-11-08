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
 * Changelog
 * ---------
 * 
 * 2014-06  Ben Dudson <benjamin.dudson@york.ac.uk>
 *     * Removed static variables in functions, changing to class members.
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
#include <utils.hxx>

#include "../../laplacefactory.hxx"

namespace {
RegisterLaplace<LaplaceSPT> registerlaplacespt(LAPLACE_SPT);
}

/// Simple parallelisation of the Thomas tridiagonal solver algorithm (serial code)
/*!
 * This is a reference code which performs the same operations as the serial code.
 * To invert a single XZ slice (FieldPerp object), data must pass from the innermost
 * processor (localmesh->PE_XIND = 0) to the outermost (localmesh->PE_XIND = localmesh->NXPE-1) and back again.
 *
 * Some parallelism is achieved by running several inversions simultaneously, so while
 * processor #1 is inverting Y=0, processor #0 is starting on Y=1. This works ok as long
 * as the number of slices to be inverted is greater than the number of X processors (MYSUB > localmesh->NXPE).
 * If MYSUB < localmesh->NXPE then not all processors can be busy at once, and so efficiency will fall sharply.
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
  LaplaceSPT(Options *opt = nullptr, const CELL_LOC = CELL_CENTRE, Mesh *mesh_in = nullptr);
  ~LaplaceSPT();
  
  using Laplacian::setCoefA;
  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Ccoef = val;
  }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceSPT does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceSPT does not have Ez coefficient");
  }

  using Laplacian::solve;
  FieldPerp solve(const FieldPerp &b) override;
  FieldPerp solve(const FieldPerp &b, const FieldPerp &x0) override;
  
  Field3D solve(const Field3D &b) override;
  Field3D solve(const Field3D &b, const Field3D &x0) override;
private:
  enum { SPT_DATA = 1123 }; ///< 'magic' number for SPT MPI messages
  
  Field2D Acoef, Ccoef, Dcoef;

  /// Data structure for SPT algorithm
  struct SPT_data {
    SPT_data() : comm_tag(SPT_DATA) {}
    void allocate(int mm, int nx); // Allocates memory
    ~SPT_data(){}; // Free memory
    
    int jy; ///< Y index
    
    Matrix<dcomplex> bk;  ///< b vector in Fourier space
    Matrix<dcomplex> xk;

    Matrix<dcomplex> gam;
  
    Matrix<dcomplex> avec, bvec, cvec; ///< Diagonal bands of matrix

    int proc; // Which processor has this reached?
    int dir;  // Which direction is it going?
  
    comm_handle recv_handle; // Handle for receives
  
    int comm_tag; // Tag for communication
  
    Array<BoutReal> buffer;
  };
  
  int ys, ye;         // Range of Y indices
  SPT_data slicedata; // Used to solve for a single FieldPerp
  SPT_data* alldata;  // Used to solve a Field3D

  Array<dcomplex> dc1d; ///< 1D in Z for taking FFTs

  void tridagForward(dcomplex *a, dcomplex *b, dcomplex *c,
                      dcomplex *r, dcomplex *u, int n,
                      dcomplex *gam,
                      dcomplex &bet, dcomplex &um, bool start=false);
  void tridagBack(dcomplex *u, int n,
                   dcomplex *gam, dcomplex &gp, dcomplex &up);
  
  int start(const FieldPerp &b, SPT_data &data);
  
  int next(SPT_data &data);
  
  void finish(SPT_data &data, FieldPerp &x);

};

#endif // __SPT_H__
