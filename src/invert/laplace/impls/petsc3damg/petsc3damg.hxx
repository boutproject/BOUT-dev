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

#ifndef __PETSC3DAMG_H__
#define __PETSC3DAMG_H__

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

class LaplacePetsc3DAmg : public Laplacian {
public:
  LaplacePetsc3DAmg(Options *opt = NULL);
  ~LaplacePetsc3DAmg(){
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

  
  const Field3D solve(const Field3D &b) override { Field3D zero(b.getMesh()); zero = 0.; return solve(b, zero); }
  const Field3D solve(const Field3D &b_in, const Field3D &x0) override;

  const FieldPerp solve(const FieldPerp& UNUSED(b)) {
    throw BoutException("LaplacePetsc3DAmg cannot solve for FieldPerp");
  }

  Field3D multiplyAx(const Field3D &x);
  
private:
  Field3D A,C1,C2,D; // ODE Coefficients
  int Nx_local, Nx_global,Nlocal,Nz_local,Nz_global,Nglobal,Ny_local,Ny_global;
  // Local and global grid sizes
  int mzstart,mxstart,mystart,lxs,lzs,lys,nxt,nzt,nyt,nxzt;
  int zbdcon,ybdcon,xbdcon;
  // 0-periodic 1-all Neumann 2-all Dirichlet
  // n == 1 mod 3: Neumann at 0
  // n == 1 mod 5: Neumann at 1

  /******* Start implementation ********/
  int mgplag,cftype,pcheck,tcheck;
  int tNP,xNP,xProcI,zNP,zProcI,yNP,yProcI;
  int xgstart,xgend,zgstart,zgend,ygstart,ygend;
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
