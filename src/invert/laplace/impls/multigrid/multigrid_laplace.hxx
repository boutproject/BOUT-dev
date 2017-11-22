/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Geometrical Multigrid Solver
 *
 * Equation solved is:
 d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2015 K.S. Kang
 * Modified version Aeptember 2015
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

#ifndef __MULTIGRID_LAPLACE_H__
#define __MULTIGRID_LAPLACE_H__

#include <mpi.h>

#include <globals.hxx>
#include <output.hxx>
#include <options.hxx>
#include <invert_laplace.hxx>
#include <boutexception.hxx>

#define MAXGM 15

// In multigrid_alg.cxx

class MultigridAlg{
public:
  MultigridAlg(int ,int ,int ,int ,int ,MPI_Comm ,int);
  virtual ~MultigridAlg();

  void setMultigridC(int );
  void getSolution(BoutReal *,BoutReal *,int ); 
  void cleanMem();

  int mglevel,mgplag,cftype,mgsm,pcheck,xNP,zNP,rProcI;
  BoutReal rtol,atol,dtol,omega;
  int *gnx,*gnz,*lnx,*lnz;
  BoutReal **matmg;

protected:
  /******* Start implementation ********/
  int numP,xProcI,zProcI,xProcP,xProcM,zProcP,zProcM;

  MPI_Comm commMG;

  void communications(BoutReal *, int );
  void setMatrixC(int );

  void cycleMG(int ,BoutReal *, BoutReal *);
  void smoothings(int , BoutReal *, BoutReal *);
  void projection(int , BoutReal *, BoutReal *);
  void prolongation(int ,BoutReal *, BoutReal *);
  void pGMRES(BoutReal *, BoutReal *, int , int);
  void solveMG(BoutReal *, BoutReal *, int );
  void multiAVec(int , BoutReal *, BoutReal *);
  void residualVec(int , BoutReal *, BoutReal *, BoutReal *);
  BoutReal vectorProd(int , BoutReal *, BoutReal *); 

  virtual void lowestSolver(BoutReal *, BoutReal *, int );
  
};


// Define three different type of multigrid solver
// in multigrid_solver.cxx

class MultigridSerial: public MultigridAlg{
public:
  MultigridSerial(int ,int ,int ,int ,int ,MPI_Comm ,int );
  ~MultigridSerial();

  void convertMatrixF(BoutReal *); 

private:
  
};

class Multigrid2DPf1D: public MultigridAlg{
public:
  Multigrid2DPf1D(int ,int ,int ,int ,int ,int ,int ,int ,MPI_Comm ,int );
  ~Multigrid2DPf1D();

  void setMultigridC(int );
  void setPcheck(int );
  void setValueS();
  void cleanS();
  int kflag;

private:
  MultigridSerial *sMG;
  void convertMatrixFS(int  ); 
  void lowestSolver(BoutReal *, BoutReal *, int );
  
};

class Multigrid1DP: public MultigridAlg{
public:
  Multigrid1DP(int ,int ,int ,int ,int ,int, MPI_Comm ,int );
  ~Multigrid1DP();
  void setMultigridC(int );
  void setPcheck(int );
  void setValueS();
  void cleanS();

  int kflag;


private:
  MPI_Comm comm2D;
  MultigridSerial *sMG;
  Multigrid2DPf1D *rMG;
  void convertMatrixF2D(int ); 
  void convertMatrixFS(int ); 
  void lowestSolver(BoutReal *, BoutReal *, int );
  
};



class LaplaceMultigrid : public Laplacian {
public:
  LaplaceMultigrid(Options *opt = NULL);
  ~LaplaceMultigrid();
  
  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C1 = val; C2 = val;  }
  void setCoefC1(const Field2D &val) { C1 = val; }
  void setCoefC2(const Field2D &val) { C2 = val; }
  void setCoefD(const Field2D &val) { D = val; }
  void setCoefEx(const Field2D &UNUSED(val)) { throw BoutException("setCoefEx is not implemented in LaplaceMultigrid"); }
  void setCoefEz(const Field2D &UNUSED(val)) { throw BoutException("setCoefEz is not implemented in LaplaceMultigrid"); }
  
  void setCoefA(const Field3D &val) { A = val; }
  void setCoefC(const Field3D &val) { C1 = val; C2 = val; }
  void setCoefC1(const Field3D &val) { C1 = val; }
  void setCoefC2(const Field3D &val) { C2 = val; }
  void setCoefD(const Field3D &val) { D = val; }
  
  const FieldPerp solve(const FieldPerp &b) { FieldPerp zero(b.getMesh()); zero = 0.; return solve(b, zero); }
  const FieldPerp solve(const FieldPerp &b_in, const FieldPerp &x0);
  
private:
  Field3D A,C1,C2,D; // ODE Coefficients
  int Nx_local, Nx_global, Nz_local, Nz_global; // Local and global grid sizes
  int yindex; // y-position of the current solution phase
  BoutReal *x; // solution vector
  BoutReal *b; // RHS vector
  Multigrid1DP *kMG;

  /******* Start implementation ********/
  int mglevel,mgplag,cftype,mgsm,pcheck,tcheck;
  int xNP,xProcI;
  int mgcount,mgmpi;

  Options *opts;
  BoutReal rtol,atol,dtol,omega;
  MPI_Comm commX;

  int comms_tagbase;

  void generateMatrixF(int);
};

#endif // __MULTIGRID_LAPLACE_H__
