/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using MUMPS Solver
 *
 **************************************************************************
 * Copyright 2013 Copyright 2013 J. Omotani (based on petsc_laplace (C) J. Buchanan)
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
#include "mumps_laplace.hxx"


#ifdef BOUT_HAS_MUMPS

// #include "mpi.h"
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <msg_stack.hxx>
#include <cmath>

LaplaceMumps::LaplaceMumps(Options *opt) : 
  Laplacian(opt),
  A(0.0), C1(1.0), C2(1.0), D(1.0), Ex(0.0), Ez(0.0),
  issetD(false), issetC(false), issetE(false)
{
  
  // Get Options in Laplace Section
  if (!opt) opts = Options::getRoot()->getSection("laplace");
  else opts=opt;
  
  #ifdef CHECK
    implemented_flags = INVERT_AC_IN_GRAD
		      + INVERT_AC_OUT_GRAD
		      + INVERT_START_NEW
// 		      + INVERT_4TH_ORDER
// 		      + INVERT_IN_SET
// 		      + INVERT_OUT_SET
		      + INVERT_IN_RHS
		      + INVERT_OUT_RHS
		      ;
    if ( flags & ~implemented_flags) {
      if (flags&INVERT_4TH_ORDER) output<<"For MUMPS based Laplacian inverter, use 'fourth_order=true' instead of setting INVERT_4TH_ORDER flag"<<endl;
      throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in petsc_laplace.cxx");
    }
  #endif

  // Get communicator for group of processors in X - all points in z-x plane for fixed y.
  comm = mesh->getXcomm();
  
  // Need to determine local size to use based on prior parallelisation
  // Coefficient values are stored only on local processors.
  localN = (mesh->xend - mesh->xstart + 1) * (mesh->ngz-1);
  if(mesh->firstX())
    localN += mesh->xstart * (mesh->ngz-1);    // If on first processor add on width of boundary region
  if(mesh->lastX())
    localN += mesh->xstart * (mesh->ngz-1);    // If on last processor add on width of boundary region
  

  // Calculate total number of points in physical grid
  if(MPI_Allreduce(&localN, &size, 1, MPI_INT, MPI_SUM, comm) != MPI_SUCCESS)
    throw BoutException("Error in MPI_Allreduce during LaplacePetsc initialisation");
  
  // Calculate total (physical) grid dimensions
  meshz = mesh->GlobalNz-1;
  meshx = size / meshz;
  
  // Calculate number of guard cells in x-direction
  nxguards = mesh->ngx - (mesh->xend-mesh->xstart+1);
  
  // Get implementation specific options
  opts->get("fourth_order", fourth_order, false);
//   opts->get("repeat_analysis", repeat_analysis, 100);

  sol.allocate();
  
  mumps_struc.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(mesh->getXcomm()); // MPI communicator for MUMPS, in fortran format
  mumps_struc.sym = 0; // Solve using unsymmetric matrix
  mumps_struc.par = 1; // Use the host processor (rank 0) to do work for the solution
  
  // Initialize the mumps struc
  mumps_struc.job = MUMPS_JOB_INIT;
  dmumps_c(&mumps_struc);
  
  // n is the order of the matrix to be inverted
  mumps_struc.n = size;
  // nz is the total number of non-zero elements in the matrix, nz_loc is the number of non-zero elements on this processor
  if (fourth_order) {
    mumps_struc.nz = 25*(meshx-nxguards)*meshz;
    mumps_struc.nz_loc = 25*(mesh->xend-mesh->xstart+1)*(mesh->ngz-1);
  }
  else {
    mumps_struc.nz = 9*(meshx-nxguards)*meshz;
    mumps_struc.nz_loc = 9*(mesh->xend-mesh->xstart+1)*(mesh->ngz-1);
  }
  if (flags & INVERT_AC_IN_GRAD) {
    if (fourth_order) {
      mumps_struc.nz += 5*meshz*mesh->xstart;
      if (mesh->firstX()) mumps_struc.nz_loc += 5*mesh->xstart*(mesh->ngz-1);
    }
    else {
      mumps_struc.nz += 3*meshz*(mesh->xstart);
      if (mesh->firstX()) mumps_struc.nz_loc += 3*mesh->xstart*(mesh->ngz-1);
    }
  }
  else {
    mumps_struc.nz += mesh->xstart*meshz;
    if (mesh->firstX()) mumps_struc.nz_loc += mesh->xstart*(mesh->ngz-1);
  }
  if (flags & INVERT_AC_OUT_GRAD) {
    if (fourth_order) {
      mumps_struc.nz += 5*(mesh->ngx-mesh->xend-1)*meshz;
      if (mesh->lastX()) mumps_struc.nz_loc += 5*(mesh->ngx-mesh->xend-1)*(mesh->ngz-1);
    }
    else {
      mumps_struc.nz += 3*(mesh->ngx-mesh->xend-1)*meshz;
      if (mesh->lastX()) mumps_struc.nz_loc += 3*(mesh->ngx-mesh->xend-1)*(mesh->ngz-1);
    }
  }
  else {
    mumps_struc.nz += (mesh->ngx-mesh->xend-1)*meshz;
    if (mesh->lastX()) mumps_struc.nz_loc += (mesh->ngx-mesh->xend-1)*(mesh->ngz-1);
  }
// // These would be needed if giving the matrix only on the host processor, or possibly if providing the structure on the host processor for analysis
//   mumps_struc.irn = new MUMPS_INT[mumps_struc.nz];
//   mumps_struc.jcn = new MUMPS_INT[mumps_struc.nz];
//   mumps_struc.a = new BoutReal[mumps_struc.nz];
  mumps_struc.irn_loc = new MUMPS_INT[mumps_struc.nz_loc]; // list of GLOBAL row indices of local matrix entries
  mumps_struc.jcn_loc = new MUMPS_INT[mumps_struc.nz_loc]; // list of GLOBAL column indices of local matrix entries
  mumps_struc.a_loc = new BoutReal[mumps_struc.nz_loc]; // the matrix entries
  
  if (mesh->firstX()) {
    mumps_struc.nrhs = 1; // number of right hand side vectors
    mumps_struc.lrhs = mumps_struc.n; // leading dimension of rhs (i.e. length of vector)
//     mumps_struc.rhs = new BoutReal[mumps_struc.lrhs*mumps_struc.nrhs]; // rhs, to be provided on the rank-0 processor only
  }
//   if (mesh->firstX()) mumps_struc.sol_loc = *sol.getData(); // pointer to the array to put the solution in, starts at 0 on first processor
//   else mumps_struc.sol_loc = *sol.getData() + mesh->xstart*meshz; // pointer to the array to put the solution in, starts at mesh->xstart
//   mumps_struc.lsol_loc = localN; // size of the (local) solution array
//   mumps_struc.isol_loc = new MUMPS_INT[localN]; // list of indices of the solution array (though this is all local points)
  mumps_struc.icntl[2] = 0; // Suppress output of global information
//   mumps_struc.icntl[3] = 0; // gives no information for debugging
  mumps_struc.icntl[3] = 1; // error messages only
//   mumps_struc.icntl[3] = 4; // gives lots of information for debugging
  opts->get("mumps_increase_working_space", mumps_struc.icntl[13], 20);
  mumps_struc.icntl[17] = 3; // provide distributed matrix for both analysis and factorization
// Don't use this: gives solution in some sructure convenient for the solver, not as desired by BOUT++       mumps_struc.icntl[20] = 1; // leave the solution distributed across the processors
//   mumps_struc.icntl[22] = ; // Apparently should provide a value here significantly larger than infog[15] when running in parallel
  mumps_struc.icntl[27] = 2; // Force parallel analysis
  mumps_struc.icntl[28] = 0; // set to 1 to force pt-scotch and to 2 to force parmetis to be used for the ordering
//   mumps_struc.cntl[0] = 0.01; // relative threshold for numerical pivoting. 0.0 is appropriate if matrix is diagonally dominant (there are no zero pivots). High values increase fill-in but lead to more accurate factorization. Default is 0.01
// for (int i=0; i<40; i++) output<<i<<" "<<mumps_struc.icntl[i]<<endl;
// exit(1);
  
//   mumps_struc.job = MUMPS_JOB_ALL;
//   iteration_count = repeat_analysis;
  
//   localrhssize = (mesh->xend-mesh->xstart+1)*mesh->ngy*mesh->ngz;
//   if (mesh->lastX()) {
//     localrhssize += (mesh->ngx-mesh->xend-1)*mesh->ngy*mesh->ngz;
//   }
//   if (mesh->firstX()) {
//     localrhssize += mesh->xstart*mesh->ngy*mesh->ngz;
//     
//     int nxpe = mesh->NXPE;
//     localrhs_size_array = new int[nxpe];
//     localrhs_size_array[0] = localrhssize;
//     if (nxpe>1) {
//       for (int i=1; i<nxpe-1; i++)
// 	localrhs_size_array[i] = (mesh->xend-mesh->xstart+1)*mesh->ngy*mesh->ngz;
//       localrhs_size_array[nxpe-1] = (mesh->ngx-mesh->xstart)*mesh->ngy*mesh->ngz;
//     }
//     rhs_positions = new int[nxpe];
//     rhs_positions[0] = 0;
//     for (int i=1; i<nxpe; i++)
//       rhs_positions[i] = rhs_positions[i-1] + localrhs_size_array[i-1];
//     
//     rhs = new BoutReal[meshx*mesh->ngy*mesh->ngz];
//     rhs_slice = new BoutReal[meshx*mesh->ngz];
//   }
  localrhssize = (mesh->xend-mesh->xstart+1)*(mesh->ngz-1);
  if (mesh->lastX()) {
    localrhssize += (mesh->ngx-mesh->xend-1)*(mesh->ngz-1);
  }
  if (mesh->firstX()) {
    localrhssize += mesh->xstart*(mesh->ngz-1);
    
    int nxpe = mesh->NXPE;
    localrhs_size_array = new int[nxpe];
    localrhs_size_array[0] = localrhssize;
    if (nxpe>1) {
      for (int i=1; i<nxpe-1; i++)
	localrhs_size_array[i] = (mesh->xend-mesh->xstart+1)*(mesh->ngz-1);
      localrhs_size_array[nxpe-1] = (mesh->ngx-mesh->xstart)*(mesh->ngz-1);
    }
    rhs_positions = new int[nxpe];
    rhs_positions[0] = 0;
    for (int i=1; i<nxpe; i++)
      rhs_positions[i] = rhs_positions[i-1] + localrhs_size_array[i-1];
    
    rhs = new BoutReal[meshx*meshz];
  }
  localrhs = new BoutReal[localrhssize];
  
  // Set Arrays of matrix indices, using i (0<=i<nz_loc), and solution indices, using j (0<=j<localN)
  int i=0; //int j=0;
  if (mesh->firstX())
    for (int x=0; x<mesh->xstart; x++)
      for (int z=0; z<mesh->ngz-1; z++) {
	int x0 = x;
	int xp = x+1;
	int xpp = x+2;
	int z0 = z;
	if(flags & INVERT_AC_IN_GRAD) {
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
	  mumps_struc.jcn_loc[i] = x0*meshz + z0 + 1;
	  i++;
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	  mumps_struc.jcn_loc[i] = xp*meshz + z0 + 1;
	  i++;
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	  mumps_struc.jcn_loc[i] = xpp*meshz + z0 + 1;
	  i++;
	  if (fourth_order) {
	    mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	    mumps_struc.jcn_loc[i] = (x+3)*meshz + z0 + 1;
	    i++;
	    mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	    mumps_struc.jcn_loc[i] = (x+4)*meshz + z0 + 1;
	    i++;
	  }
	}
	else {
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
	  mumps_struc.jcn_loc[i] = x0*meshz + z0 + 1;
	  i++;
	}
// 	mumps_struc.isol_loc[j] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
// 	j++;
      }
  for (int x=mesh->xstart; x<=mesh->xend; x++)
    for (int z=0; z<mesh->ngz-1; z++) {
      int xmm = mesh->XGLOBAL(x)-2;
      int xm = mesh->XGLOBAL(x)-1;
      int x0 = mesh->XGLOBAL(x);
      int xp = mesh->XGLOBAL(x)+1;
      int xpp = mesh->XGLOBAL(x)+2;
      int zmm = (z-2<0) ? (z-2+meshz) : (z-2);
      int zm = (z-1<0) ? (z-1+meshz) : (z-1);
      int z0 = z;
      int zp = (z+1>=meshz) ? (z+1-meshz) : (z+1);
      int zpp = (z+2>=meshz) ? (z+2-meshz) : (z+2);
      if (fourth_order) {
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
	mumps_struc.jcn_loc[i] = x0*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xmm*meshz + zmm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xmm*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xmm*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xmm*meshz + (z+1) + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xmm*meshz + zpp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + zmm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + zp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + zpp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = x0*meshz + zmm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = x0*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = x0*meshz + zp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = x0*meshz + zpp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + zmm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + zp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + zpp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xpp*meshz + zmm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xpp*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xpp*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xpp*meshz + zp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xpp*meshz + zpp + 1;
	i++;
      }
      else {
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
	mumps_struc.jcn_loc[i] = x0*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xm*meshz + zp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = x0*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = x0*meshz + zp + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + zm + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + z0 + 1;
	i++;
	mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	mumps_struc.jcn_loc[i] = xp*meshz + zp + 1;
	i++;
      }
//       mumps_struc.isol_loc[j] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
//       j++;
    }
  if (mesh->lastX())
    for (int x=mesh->xend+1; x<mesh->ngx; x++)
      for (int z=0; z<mesh->ngz-1; z++) {
	int xmm = mesh->XGLOBAL(mesh->xend)+x-mesh->xend-2;
	int xm = mesh->XGLOBAL(mesh->xend)+x-mesh->xend-1;
	int x0 = mesh->XGLOBAL(mesh->xend)+x-mesh->xend;
	int z0 = z;
	if(flags & INVERT_AC_IN_GRAD) {
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
	  mumps_struc.jcn_loc[i] = x0*meshz + z0 + 1;
	  i++;
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	  mumps_struc.jcn_loc[i] = xm*meshz + z0 + 1;
	  i++;
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	  mumps_struc.jcn_loc[i] = xmm*meshz + z0 + 1;
	  i++;
	  if (fourth_order) {
	    mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	    mumps_struc.jcn_loc[i] = (x0-3)*meshz + z0 + 1;
	    i++;
	    mumps_struc.irn_loc[i] = x0*meshz + z0 + 1;
	    mumps_struc.jcn_loc[i] = (x0-4)*meshz + z0 + 1;
	    i++;
	  }
	}
	else {
	  mumps_struc.irn_loc[i] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
	  mumps_struc.jcn_loc[i] = x0*meshz + z0 + 1;
	  i++;
	}
// 	mumps_struc.isol_loc[j] = x0*meshz + z0 + 1; // Indices for fortran arrays that start at 1
// 	j++;
      }
  
  if ( i!=mumps_struc.nz_loc ) bout_error("LaplaceMumps: matrix index error");
//   if ( j!=localN ) bout_error("LaplaceMumps: vector index error");
// output<<"matrix indices:"<<endl;for (int k=0; k<mumps_struc.nz; k++) output<<k<<" "<<mumps_struc.irn_loc[k]<<" "<<mumps_struc.jcn_loc[k]<<endl;
// output<<endl<<"solution vector indices:"<<endl;for (int k=0; k<mumps_struc.n;k++) output<<k<<" "<<mumps_struc.isol_loc[k]<<endl;
// output<<"nz="<<mumps_struc.nz<<" nz_loc="<<mumps_struc.nz_loc<<endl;
// MPI_Barrier(BoutComm::get()); exit(13);
  
  mumps_struc.job = MUMPS_JOB_ANALYSIS;
  dmumps_c( &mumps_struc );
  mumps_struc.job = MUMPS_JOB_BOTH;
}

// const Field3D LaplaceMumps::solve(const Field3D &b, const Field3D &x0) {
//   return solve(b);
// }
// 
// const Field3D LaplaceMumps::solve(const Field3D &b) {
//   Timer timer("invert");
// #ifdef CHECK
//   msg_stack.push("Laplacian::solve(Field3D)");
// #endif
//   int ys = mesh->ystart, ye = mesh->yend;
// 
//   if(mesh->hasBndryLowerY()) {
//     if (include_yguards)
//       ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells
//     
//     ys += extra_yguards_lower;
//   }
//   if(mesh->hasBndryUpperY()) {
//     if (include_yguards)
//       ye = mesh->ngy-1; // Contains upper boundary and we are solving in the guard cells
//       
//     ye -= extra_yguards_upper;
//   }
// 
//   Field3D x = copy(b); // Force new memory allocation as we will mess around with x's data via pointers (i.e. 'unsafely')
//   
//   BoutReal* localrhs = **x.getData(); // Input the rhs in the solution field as solution will be returned in place by MUMPS
//   if (!mesh->firstX()) localrhs += mesh->xstart*mesh->ngy*mesh->ngz;
//   MPI_Gatherv(localrhs,localrhssize,MPI_DOUBLE,rhs,localrhs_size_array,rhs_positions,MPI_DOUBLE,0,mesh->getXcomm());
//   
//   if ( ++iteration_count > repeat_analysis ) {
//     mumps_struc.job = MUMPS_JOB_ALL;
//     for(int jy=ys; jy <= ye; jy++) {
//       if (mesh->firstX())
// 	for(int jx=0; jx<meshx; jx++)
// 	  for (int jz=0; jz<meshz; jz++)
// 	    rhs_slice[jx*meshz+jz] = rhs[jx*mesh->ngy*mesh->ngz + jy*mesh->ngz + jz];
//       
//       solve(rhs_slice,jy);
//       
//       if (mesh->firstX())
// 	for(int jx=0; jx<meshx; jx++)
// 	  for (int jz=0; jz<meshz; jz++)
// 	     rhs[jx*mesh->ngy*mesh->ngz + jy*mesh->ngz + jz] = rhs_slice[jx*meshz+jz];
//     }
//     mumps_struc.job = MUMPS_JOB_BOTH;
//   }
//   else {
//     for(int jy=ys; jy <= ye; jy++) {
//       if (mesh->firstX())
// 	for(int jx=0; jx<meshx; jx++)
// 	  for (int jz=0; jz<meshz; jz++)
// 	    rhs_slice[jx*meshz+jz] = rhs[jx*mesh->ngy*mesh->ngz + jy*mesh->ngz + jz];
//       
//       solve(rhs_slice,jy);
//       
//       if (mesh->firstX())
// 	for(int jx=0; jx<meshx; jx++)
// 	  for (int jz=0; jz<meshz; jz++)
// 	     rhs[jx*mesh->ngy*mesh->ngz + jy*mesh->ngz + jz] = rhs_slice[jx*meshz+jz];
//     }
//   }
//   
//   MPI_Scatterv(rhs,localrhs_size_array,rhs_positions,MPI_DOUBLE,localrhs,localrhssize,MPI_DOUBLE,0,mesh->getXcomm()); // Scatters solution from host back to localrhs (which points to x's data) on all processors
//   
// #ifdef CHECK
//   msg_stack.pop();
// #endif
// 
//   x.setLocation(b.getLocation());
// 
//   return x;
// }

// const Field3D LaplaceMumps::solve(const Field3D &b, const Field3D &x0) {
//   return solve(b);
// }
// 
// const Field3D LaplaceMumps::solve(const Field3D &b) {
//   Timer timer("invert");
// #ifdef CHECK
//   msg_stack.push("Laplacian::solve(Field3D)");
// #endif
//   int ys = mesh->ystart, ye = mesh->yend;
// 
//   if(mesh->hasBndryLowerY()) {
//     if (include_yguards)
//       ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells
//     
//     ys += extra_yguards_lower;
//   }
//   if(mesh->hasBndryUpperY()) {
//     if (include_yguards)
//       ye = mesh->ngy-1; // Contains upper boundary and we are solving in the guard cells
//       
//     ye -= extra_yguards_upper;
//   }
//   
//   Field3D x;
//   x.allocate();
//   
//   if ( ++iteration_count > repeat_analysis ) {
//     mumps_struc.job = MUMPS_JOB_ALL;
//     for(int jy=ys; jy <= ye; jy++) {
//       x = solve(b.slice(jy));
//     }
//     mumps_struc.job = MUMPS_JOB_BOTH;
//     iteration_count = 0;
//   }
//   else {
//     for(int jy=ys; jy <= ye; jy++) {
//       x = solve(b.slice(jy));
//     }
//   }
//   
// #ifdef CHECK
//   msg_stack.pop();
// #endif
// 
//   x.setLocation(b.getLocation());
// 
//   return x;
// 
// }

const FieldPerp LaplaceMumps::solve(const FieldPerp &b, const FieldPerp &x0) {
  return solve(b);
}

const FieldPerp LaplaceMumps::solve(const FieldPerp &b) {
  int y = b.getIndex();
  sol = 0.;
  sol.setIndex(y);
  
  // Set Dirichlet boundary conditions through rhs if needed
  if (!(flags && INVERT_AC_IN_GRAD+INVERT_IN_RHS)) {
    if (mesh->firstX())
      for (int z=0; z<mesh->ngz-1; z++)
	for (int x=mesh->xstart-1; x>=0; x--) {
	  b[x][z]=0.;
	}
  }
  if (!(flags && INVERT_AC_OUT_GRAD+INVERT_OUT_RHS)) {
    if (mesh->lastX())
      for (int z=0; z<mesh->ngz-1; z++)
	for (int x=mesh->xend+1; x<mesh->ngx; x++) {
	  b[x][z]=0.;
	}
  }
  
  BoutReal* bdata = *b.getData();
  int xs,xe;
  if (mesh->firstX()) xs=0;
  else xs=mesh->xstart;
  if (mesh->lastX()) xe=mesh->ngx-1;
  else xe=mesh->xend;
  
  for (int x=xs; x<=xe; x++)
    for (int z=0; z<mesh->ngz-1; z++)
      localrhs[(x-xs)*(mesh->ngz-1)+z] = bdata[x*mesh->ngz+z];
  
  MPI_Gatherv(localrhs,localrhssize,MPI_DOUBLE,rhs,localrhs_size_array,rhs_positions,MPI_DOUBLE,0,comm);
  
  solve(rhs,y);
  
  MPI_Scatterv(rhs,localrhs_size_array,rhs_positions,MPI_DOUBLE,localrhs,localrhssize,MPI_DOUBLE,0,comm); // Scatters solution from host back to localrhs (which points to x's data) on all processors
  
  BoutReal* soldata = *sol.getData();
  for (int x=xs; x<=xe; x++)
    for (int z=0; z<mesh->ngz-1; z++)
      soldata[x*mesh->ngz+z] = localrhs[(x-xs)*(mesh->ngz-1)+z];
  
  return sol;
}

void LaplaceMumps::solve(BoutReal* rhs, int y) {

{ Timer timer("mumpssetup");
  int i = 0;
  
  // Set Matrix Elements corresponding to index lists created in constructor (x,z) loop over rows

  // X=0 to mesh->xstart-1 defines the boundary region of the domain.
  if( mesh->firstX() )
    for(int x=0; x<mesh->xstart; x++)
      for(int z=0; z<mesh->ngz-1; z++) {
	// Set values corresponding to nodes adjacent in x if Neumann Boundary Conditions are required.
	if(flags & INVERT_AC_IN_GRAD)
	  if( fourth_order ) {
	    // Fourth Order Accuracy on Boundary
	    mumps_struc.a_loc[i] = -25.0 / (12.0*mesh->dx[x][y]);
	    i++;
	    mumps_struc.a_loc[i] = 4.0 / mesh->dx[x][y];
	    i++;
	    mumps_struc.a_loc[i] = -3.0 / mesh->dx[x][y];
	    i++;
	    mumps_struc.a_loc[i] = 4.0 / (3.0*mesh->dx[x][y]);
	    i++;
	    mumps_struc.a_loc[i] = -1.0 / (4.0*mesh->dx[x][y]);
	    i++;
	  }
	  else {
	    // Second Order Accuracy on Boundary
	    mumps_struc.a_loc[i] = -3.0 / (2.0*mesh->dx[x][y]);
	    i++;
	    mumps_struc.a_loc[i] = 2.0 / mesh->dx[x][y];
	    i++;
	    mumps_struc.a_loc[i] = -1.0 / (2.0*mesh->dx[x][y]);
	    i++;
	  }
	else {
	  // Set Diagonal Values to 1
	  mumps_struc.a_loc[i] = 1.;
	  i++;
	}
      }
  
  // Main domain with Laplacian operator
  for(int x=mesh->xstart; x <= mesh->xend; x++)
    for(int z=0; z<mesh->ngz-1; z++) {
      BoutReal A0, A1, A2, A3, A4, A5;
      A0 = A[x][y][z];
      Coeffs( x, y, z, A1, A2, A3, A4, A5 );
      
      BoutReal dx   = mesh->dx[x][y];
      BoutReal dx2  = pow( mesh->dx[x][y] , 2.0 );
      BoutReal dz   = mesh->dz;
      BoutReal dz2  = pow( mesh->dz, 2.0 );
      BoutReal dxdz = mesh->dx[x][y] * mesh->dz;
      
      // Set Matrix Elements
      if (fourth_order) {
	// f(i,j) = f(x,z)
	mumps_struc.a_loc[i] = A0 - (5.0/2.0)*( (A1 / dx2) + (A2 / dz2) );
	i++;
	
	// f(i-2,j-2)
	mumps_struc.a_loc[i] = A3 / ( 144.0 * dxdz );
	i++;

	// f(i-2,j-1)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 18.0 * dxdz );
	i++;

	// f(i-2,j)
	mumps_struc.a_loc[i] = (1.0/12.0) * ( (-1.0 * A1 /  dx2 ) + (A4 / dx) ); 
	i++;

	// f(i-2,j+1)
	mumps_struc.a_loc[i] = A3 / ( 18.0 * dxdz );
	i++;

	// f(i-2,j+2)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 144.0 * dxdz );
	i++;

	// f(i-1,j-2)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 18.0 * dxdz );
	i++;

	// f(i-1,j-1)
	mumps_struc.a_loc[i] = 4.0 * A3 / ( 9.0 * dxdz ); 
	i++;

	// f(i-1,j)
	mumps_struc.a_loc[i] = ( 4.0 * A1 / ( 3.0 * dx2 ) ) - ( 2.0 * A4 / ( 3.0 * dx ) );
	i++;

	// f(i-1,j+1)
	mumps_struc.a_loc[i] = -4.0 * A3 / ( 9.0 * dxdz ); 
	i++;

	// f(i-1,j+2)
	mumps_struc.a_loc[i] = A3 / ( 18.0 * dxdz );
	i++;

	// f(i,j-2)
	mumps_struc.a_loc[i] = (1.0/12.0) * ( ( -1.0 * A2 / dz2 ) + ( A5 / dz ) ); 
	i++;
	
	// f(i,j-1)
	mumps_struc.a_loc[i] = ( 4.0 * A2 / ( 3.0 * dz2 ) ) - ( 2.0 * A5 / ( 3.0 * dz ) );
	i++;

	// f(i,j+1)
	mumps_struc.a_loc[i] = ( 4.0 * A2 / ( 3.0 * dz2 ) ) + ( 2.0 * A5 / ( 3.0 * dz ) );
	i++;

	// f(i,j+2)
	mumps_struc.a_loc[i] = (-1.0/12.0) * ( ( A2 / dz2 ) + ( A5 / dz ) ); 
	i++;

	// f(i+1,j-2)
	mumps_struc.a_loc[i] = A3 / ( 18.0 * dxdz );
	i++;
	
	// f(i+1,j-1)
	mumps_struc.a_loc[i] = -4.0 * A3 / ( 9.0 * dxdz );
	i++;
	
	// f(i+1,j)
	mumps_struc.a_loc[i] = ( 4.0 * A1 / ( 3.0*dx2 ) ) + ( 2.0 * A4 / ( 3.0 * dx ) ); 
	i++;
	
	// f(i+1,j+1)
	mumps_struc.a_loc[i] = 4.0 * A3 / ( 9.0 * dxdz ); 
	i++;

	// f(i+1,j+2)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 18.0 * dxdz );
	i++;

	// f(i+2,j-2)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 144.0 * dxdz );
	i++;

	// f(i+2,j-1)
	mumps_struc.a_loc[i] = A3 / ( 18.0 * dxdz );
	i++;

	// f(i+2,j)
	mumps_struc.a_loc[i] = (-1.0/12.0) * ( (A1 / dx2) + (A4 / dx) ); 
	i++;

	// f(i+2,j+1)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 18.0 * dxdz );
	i++;

	// f(i+2,j+2)
	mumps_struc.a_loc[i] = A3 / ( 144.0 * dxdz );
	i++;
      }
      else {
	// f(i,j) = f(x,z)
	mumps_struc.a_loc[i] = A0 - 2.0*( (A1 / dx2) + (A2 / dz2) ); 
	i++;
	
	// f(i-1,j-1)
	mumps_struc.a_loc[i] = A3 / (4.0 * dxdz); 
	i++;

	// f(i-1,j)
	mumps_struc.a_loc[i] = ( A1 / dx2 ) - A4 / ( 2.0 * dx ); 
	i++;

	// f(i-1,j+1)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 4.0 * dxdz ); 
	i++;

	// f(i,j-1)
	mumps_struc.a_loc[i] = ( A2 / dz2 ) - ( A5 / ( 2.0 * dz ) ); 
	i++;

	// f(i,j+1)
	mumps_struc.a_loc[i] = ( A2 / dz2 ) + ( A5 / ( 2.0 * dz ) ); 
	i++;

	// f(i+1,j-1)
	mumps_struc.a_loc[i] = -1.0 * A3 / ( 4.0 * dxdz ); 
	i++;
	
	// f(i+1,j)
	mumps_struc.a_loc[i] = ( A1 / dx2 ) + ( A4 / ( 2.0 * dx ) ); 
	i++;
	
	// f(i+1,j+1)
	mumps_struc.a_loc[i] = A3 / ( 4.0 * dxdz ); 
	i++;
      }
    }
  
  // X=mesh->xend+1 to mesh->ngx-1 defines the upper boundary region of the domain.
  if( mesh->lastX() )
    for(int x=mesh->xend+1; x<mesh->ngx; x++) 
      for(int z=0; z<mesh->ngz-1; z++) {
	
	// Set values corresponding to nodes adjacent in x if Neumann Boundary Conditions are required.
	if(flags & INVERT_AC_OUT_GRAD) {
	  if( fourth_order ) {
	    // Fourth Order Accuracy on Boundary
	    mumps_struc.a_loc[i] = 25.0 / (12.0*mesh->dx[x][y]); 
	    i++;
	    mumps_struc.a_loc[i] = -4.0 / mesh->dx[x][y]; 
	    i++;
	    mumps_struc.a_loc[i] = 3.0 / mesh->dx[x][y]; 
	    i++;
	    mumps_struc.a_loc[i] = -4.0 / (3.0*mesh->dx[x][y]); 
	    i++;
	    mumps_struc.a_loc[i] = 1.0 / (4.0*mesh->dx[x][y]); 
	    i++;
	  }
	  else {
	    // Second Order Accuracy on Boundary
	    mumps_struc.a_loc[i] = 3.0 / (2.0*mesh->dx[x][y]); 
	    i++;
	    mumps_struc.a_loc[i] = -2.0 / mesh->dx[x][y]; 
	    i++;
	    mumps_struc.a_loc[i] = 1.0 / (2.0*mesh->dx[x][y]); 
	    i++;
	  }
	}
	else {
	  // Set Diagonal Values to 1
	  mumps_struc.a_loc[i] = 1; 
	  i++;
	}
      }
  
  if ( i!=mumps_struc.nz_loc ) bout_error("LaplaceMumps: matrix index error");
  
  mumps_struc.rhs = rhs;
}
{ Timer timer("mumpssolve");
  // Solve the system
  dmumps_c( &mumps_struc );
}

}

void LaplaceMumps::Coeffs( int x, int y, int z, BoutReal &coef1, BoutReal &coef2, BoutReal &coef3, BoutReal &coef4, BoutReal &coef5 )
{
  coef1 = mesh->g11[x][y];     // X 2nd derivative coefficient
  coef2 = mesh->g33[x][y];     // Z 2nd derivative coefficient
  coef3 = 2.*mesh->g13[x][y];  // X-Z mixed derivative coefficient
  
  coef4 = 0.0;
  coef5 = 0.0;
  if(all_terms) {
    coef4 = mesh->G1[x][y]; // X 1st derivative
    coef5 = mesh->G3[x][y]; // Z 1st derivative
  }
  
  if(nonuniform) 
    {
      // non-uniform mesh correction
      if((x != 0) && (x != (mesh->ngx-1))) 
	{
	  //coef4 += mesh->g11[jx][jy]*0.25*( (1.0/dx[jx+1][jy]) - (1.0/dx[jx-1][jy]) )/dx[jx][jy]; // SHOULD BE THIS (?)
	  //coef4 -= 0.5 * ( ( mesh->dx[x+1][y] - mesh->dx[x-1][y] ) / SQ ( mesh->dx[x][y] ) ) * coef1; // BOUT-06 term
	  coef4 -= 0.5 * ( ( mesh->dx[x+1][y] - mesh->dx[x-1][y] ) / pow( mesh->dx[x][y], 2.0 ) ) * coef1; // BOUT-06 term
	}
    }
  
  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // d2dz2 term
    coef2 += mesh->g11[x][y] * mesh->IntShiftTorsion[x][y] * mesh->IntShiftTorsion[x][y];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }
  
  if (issetD) {
    coef1 *= D[x][y][z];
    coef2 *= D[x][y][z];
    coef3 *= D[x][y][z];
    coef4 *= D[x][y][z];
    coef5 *= D[x][y][z];
  }
  
  // A second/fourth order derivative term
  if (issetC) {
//     if( (x > 0) && (x < (mesh->ngx-1)) ) //Valid if doing second order derivative, not if fourth: should only be called for xstart<=x<=xend anyway
    if( (x > 1) && (x < (mesh->ngx-2)) ) {
      int zp = z+1;
      if (zp > meshz-1) zp -= meshz;
      int zm = z-1;
      if (zm<0) zm += meshz;
      BoutReal ddx_C;
      BoutReal ddz_C;
      if (fourth_order) {
	int zpp = z+2;
	if (zpp > meshz-1) zpp -= meshz;
	int zmm = z-2;
	if (zmm<0) zmm += meshz;
	ddx_C = (-C2[x+2][y][z] + 8.*C2[x+1][y][z] - 8.*C2[x-1][y][z] + C2[x-2][y][z]) / (12.*mesh->dx[x][y]*(C1[x][y][z]));
	ddz_C = (-C2[x][y][zpp] + 8.*C2[x][y][zp] - 8.*C2[x][y][zm] + C2[x][y][zmm]) / (12.*mesh->dz*(C1[x][y][z]));
      }
      else {
	ddx_C = (C2[x+1][y][z] - C2[x-1][y][z]) / (2.*mesh->dx[x][y]*(C1[x][y][z]));
	ddz_C = (C2[x][y][zp] - C2[x][y][zm]) / (2.*mesh->dz*(C1[x][y][z]));
      }
      
      coef4 += mesh->g11[x][y] * ddx_C + mesh->g13[x][y] * ddz_C;
      coef5 += mesh->g13[x][y] * ddx_C + mesh->g33[x][y] * ddz_C;
    }
  }
  
  // Additional 1st derivative terms to allow for solution field to be component of vector
  // NB multiply by D or Grad_perp(C)/C as appropriate before passing to setCoefEx()/setCoefEz() because both (in principle) are needed and we don't know how to split them up here
  if (issetE) {
    coef4 += Ex[x][y][z];
    coef5 += Ez[x][y][z];
  }
  
}

#endif // BOUT_HAS_MUMPS
