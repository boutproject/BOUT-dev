/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Geometrical Multigrid Solver
 *
 * Equation solved is:
 d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2018 K.S. Kang
 * Modified version September 2018
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

#ifndef __PETSCAMG_H__
#define __PETSCAMG_H__

#include <mpi.h>

#include <globals.hxx>
#include <output.hxx>
#include <options.hxx>
#include <invert_laplace.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <bout/openmpwrap.hxx>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef BOUT_HAS_PETSC

#else

// #include <petscksp.h>
#include <petscpc.h>
#include <bout/petsclib.hxx>
#include <msg_stack.hxx>
#include <bout/openmpwrap.hxx>

#ifdef _OPENMP
#include <omp.h>
#endif

class LaplacePetscAmg : public Laplacian {
public:
  LaplacePetscAmg(Options *opt = NULL);
  ~LaplacePetscAmg(){
    ISLocalToGlobalMappingDestroy(&mgmapping);
    VecDestroy( &xs );
    VecDestroy( &bs );
    //    delete [] ksp;
    //    delete [] pc;
    delete [] gindices;
  }
  
  void setCoefA(const Field2D &val) override { A = val; }
  void setCoefC(const Field2D &val) override { C1 = val; C2 = val;  }
  void setCoefC1(const Field2D &val) override { C1 = val; }
  void setCoefC2(const Field2D &val) override { C2 = val; }
  void setCoefD(const Field2D &val) override { D = val; }
  void setCoefEx(const Field2D &UNUSED(val)) override { throw BoutException("setCoefEx is not implemented in LaplaceMultigrid"); }
  void setCoefEz(const Field2D &UNUSED(val)) override { throw BoutException("setCoefEz is not implemented in LaplaceMultigrid"); }
  
  void setCoefA(const Field3D &val) override { A = val; }
  void setCoefC(const Field3D &val) override { C1 = val; C2 = val; }
  void setCoefC1(const Field3D &val) override { C1 = val; }
  void setCoefC2(const Field3D &val) override { C2 = val; }
  void setCoefD(const Field3D &val) override { D = val; }
  void settingSolver(int);

  
  FieldPerp solve(const FieldPerp &b) override { FieldPerp zero(b.getMesh()); zero = 0.; return solve(b, zero); }
  FieldPerp solve(const FieldPerp &b_in, const FieldPerp &x0) override;

  FieldPerp multiplyAx(const FieldPerp &x);
  
private:
  Field3D A,C1,C2,D; // ODE Coefficients
  int Nx_local, Nx_global,Nlocal,Nz_local,Nz_global,Nglobal;
  // Local and global grid sizes
  int mzstart,mxstart,lxs,lzs,nxt,nzt;
  int yindex; // y-position of the current solution phase

  /******* Start implementation ********/
  int mgplag,cftype,pcheck,tcheck;
  int xNP,xProcI,zNP,zProcI,xgstart,xgend,zgstart,zgend;
  int mgcount,mgmpi;

  /* Options for solver */
  Options *opts;
  BoutReal rtol,atol,dtol,omega;
  bool rightpre;
  int diffpre,elemf;
  int maxits,mgsm,mglevel,fcheck,stypei;
  std::string soltype;  /* direct, gmres, gamg, boomermg, ml, kmg, ...  */

  /*********************************************************/
  PetscLib lib;
  
  MPI_Comm commX;
  Mat MatA,MatP;
  Vec xs, bs;                 // Solution and RHS vectors
  KSP ksp;
  PC pc;
  std::string ksptype, pctype;	      // Preconditioner type
  ISLocalToGlobalMapping mgmapping;
  
  int comms_tagbase,*gindices;

  void generateMatrixA(int);
  void generateMatrixP(int);
};
#endif //   BOUT_HAS_PETSC
#endif // __PETSCAMG_H__
