/**************************************************************************
 * Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is: d*\nabla^2_\perp x + (1/c)\nabla_perp c\cdot\nabla_\perp x + a x = b
 * 
 * Where a, c and d are functions of x and y only (not z)
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

class Laplacian;

#ifndef __LAPLACE_H__
#define __LAPLACE_H__

#ifdef BOUT_HAS_PETSC
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#endif

#include "fieldperp.hxx"
#include "field3d.hxx"
#include "field2d.hxx"

#include "dcomplex.hxx"
#include "options.hxx"

// Inversion flags

const int INVERT_DC_IN_GRAD  = 1;
const int INVERT_AC_IN_GRAD  = 2;
const int INVERT_DC_OUT_GRAD = 4;
const int INVERT_AC_OUT_GRAD = 8;
const int INVERT_ZERO_DC     = 16;
const int INVERT_START_NEW   = 32;
const int INVERT_BNDRY_ONE   = 64; // Sets the width of the boundary to 1
const int INVERT_4TH_ORDER   = 128; // Use band solver for 4th order in x

const int INVERT_AC_IN_LAP   = 256;
const int INVERT_AC_OUT_LAP  = 512;

const int INVERT_IN_SYM  =  1024; // Use symmetry to enforce either zero-value or zero-gradient
const int INVERT_OUT_SYM =  2048; // Same for outer boundary
const int INVERT_IN_SET  =  4096; // Set inner boundary
const int INVERT_OUT_SET =  8192; // Set outer boundary
const int INVERT_IN_RHS  = 16384; // Use input value in RHS at inner boundary
const int INVERT_OUT_RHS = 32768; // Use input value in RHS at outer boundary
const int INVERT_KX_ZERO = 65536; // Zero the kx=0, n = 0 component

const int INVERT_DC_IN_LAP = 131072;

const int INVERT_BNDRY_IN_ONE = 262144;
const int INVERT_BNDRY_OUT_ONE = 524288;
const int INVERT_DC_IN_GRADPAR = 1048576;
const int INVERT_DC_IN_GRADPARINV = 2097152;

/// Base class for Laplacian inversion
class Laplacian {
public:
  Laplacian(Options *options = NULL);
  
  /// Set coefficients for inversion. Re-builds matrices if necessary
  virtual void setCoefA(const Field2D &val) = 0;
  virtual void setCoefA(const Field3D &val) { setCoefA(val.DC()); }
  virtual void setCoefA(const BoutReal &r) { Field2D f(r); setCoefA(f); }
  
  virtual void setCoefC(const Field2D &val) = 0;
  virtual void setCoefC(const Field3D &val) { setCoefC(val.DC()); }
  virtual void setCoefC(const BoutReal &r) { Field2D f(r); setCoefC(f); }
  
  virtual void setCoefD(const Field2D &val) = 0;
  virtual void setCoefD(const Field3D &val) { setCoefD(val.DC()); }
  virtual void setCoefD(const BoutReal &r) { Field2D f(r); setCoefD(f); }
  
  virtual void setCoefEx(const Field2D &val) = 0;
  virtual void setCoefEx(const Field3D &val) { setCoefEx(val.DC()); }
  virtual void setCoefEx(const BoutReal &r) { Field2D f(r); setCoefEx(f); }
  
  virtual void setCoefEz(const Field2D &val) = 0;
  virtual void setCoefEz(const Field3D &val) { setCoefEz(val.DC()); }
  virtual void setCoefEz(const BoutReal &r) { Field2D f(r); setCoefD(f); }
  
  virtual void setFlags(int f) {flags = f;}
  
  virtual const FieldPerp solve(const FieldPerp &b) = 0;
  virtual const Field3D solve(const Field3D &b);
  virtual const Field2D solve(const Field2D &b);
  
  virtual const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0) { return solve(b); }
  virtual const Field3D solve(const Field3D &b, const Field3D &x0);
  virtual const Field2D solve(const Field2D &b, const Field2D &x0);

  /// Coefficients in tridiagonal inversion
  void tridagCoefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c, const Field2D *ccoef = NULL, const Field2D *d=NULL);
  
  static Laplacian* create(Options *opt = NULL);  ///< Create a new Laplacian solver
  static Laplacian* defaultInstance(); ///< Return pointer to global singleton
  
protected:
  bool async_send; ///< If true, use asyncronous send in parallel algorithms
  
  int maxmode;     ///< The maximum Z mode to solve for
  
  bool low_mem;    ///< If true, reduce the amount of memory used
  bool all_terms;  // applies to Delp2 operator and laplacian inversion
  bool nonuniform; // Non-uniform mesh correction
  bool include_yguards; // solve in y-guard cells, default true.
  int extra_yguards_lower; // exclude some number of points at the lower boundary, useful for staggered grids or when boundary conditions make inversion redundant
  int extra_yguards_upper; // exclude some number of points at the upper boundary, useful for staggered grids or when boundary conditions make inversion redundant
  
  int flags;       ///< Default flags

  void tridagMatrix(dcomplex **avec, dcomplex **bvec, dcomplex **cvec,
                    dcomplex **bk, int jy, int flags, 
                    const Field2D *a = NULL, const Field2D *ccoef=NULL, 
                    const Field2D *d = NULL);
private:
  /// Singleton instance
  static Laplacian *instance;
};

////////////////////////////////////////////
// Legacy interface
// These will be removed at some point

void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c, const Field2D *ccoef = NULL, const Field2D *d=NULL);

int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c=NULL, const Field2D *d=NULL);
int invert_laplace(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c=NULL, const Field2D *d=NULL);

/// More readable API for calling Laplacian inversion. Returns x
const Field3D invert_laplace(const Field3D &b, int flags, 
                             const Field2D *a = NULL, const Field2D *c=NULL, const Field2D *d=NULL);


#endif // __LAPLACE_H__

