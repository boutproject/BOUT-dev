/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using PETSc Solvers
 *
 **************************************************************************
 * Copyright 2013 J. Buchanan, J.Omotani
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
#ifdef BOUT_HAS_PETSC_3_3

#include "petsc_laplace.hxx"

#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>

#define KSP_RICHARDSON  "richardson"
#define KSP_CHEBYSHEV   "chebyshev"
#define KSP_CG          "cg"
#define KSP_GMRES       "gmres"
#define KSP_TCQMR       "tcqmr"
#define KSP_BCGS        "bcgs"
#define KSP_CGS         "cgs"
#define KSP_TFQMR       "tfqmr"
#define KSP_CR          "cr"
#define KSP_LSQR        "lsqr"
#define KSP_BICG        "bicg"
#define KSP_PREONLY     "preonly"

LaplacePetsc::LaplacePetsc(Options *opt) : 
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
		      + INVERT_IN_SET
		      + INVERT_OUT_SET
		      + INVERT_IN_RHS
		      + INVERT_OUT_RHS
		      ;
    if ( flags & ~implemented_flags) {
      if (flags&INVERT_4TH_ORDER) output<<"For PETSc based Laplacian inverter, use 'fourth_order=true' instead of setting INVERT_4TH_ORDER flag"<<endl;
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
  meshz = mesh->ngz-1;
  meshx = size / meshz;

  // Create Vectors 
  VecCreate( comm, &xs );                        
  VecSetSizes( xs, localN, size );                
  VecSetFromOptions( xs );                     
  VecDuplicate( xs , &bs );                   
  
  // Get 4th order solver switch
  opts->get("fourth_order", fourth_order, false);
  
  // Set size of Matrix on each processor to localN x localN
  MatCreate( comm, &MatA );                                
  MatSetSizes( MatA, localN, localN, size, size );
  MatSetFromOptions(MatA);
  //   if (fourth_order) MatMPIAIJSetPreallocation( MatA, 25, PETSC_NULL, 10, PETSC_NULL );
//   else MatMPIAIJSetPreallocation( MatA, 9, PETSC_NULL, 3, PETSC_NULL );
  PetscInt *d_nnz, *o_nnz;
  PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
  PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );
  if (fourth_order) {
    // first and last 2*(mesh-ngz-1) entries are the edge x-values that (may) have 'off-diagonal' components (i.e. on another processor)
    if ( mesh->firstX() && mesh->lastX() ) {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=15;
	d_nnz[localN-1-i]=15;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=0;
      }
      for (int i=(mesh->ngz-1); i<2*(mesh->ngz-1); i++) {
	d_nnz[i]=20;
	d_nnz[localN-1-i]=20;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=0;
      }
    }
    else if ( mesh->firstX() ) {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=15;
	d_nnz[localN-1-i]=15;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=10;
      }
      for (int i=(mesh->ngz-1); i<2*(mesh->ngz-1); i++) {
	d_nnz[i]=20;
	d_nnz[localN-1-i]=20;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=5;
      }
    }
    else if ( mesh->lastX() ) {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=15;
	d_nnz[localN-1-i]=15;
	o_nnz[i]=10;
	o_nnz[localN-1-i]=0;
      }
      for (int i=(mesh->ngz-1); i<2*(mesh->ngz-1); i++) {
	d_nnz[i]=20;
	d_nnz[localN-1-i]=20;
	o_nnz[i]=5;
	o_nnz[localN-1-i]=0;
      }
    }
    else {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=15;
	d_nnz[localN-1-i]=15;
	o_nnz[i]=10;
	o_nnz[localN-1-i]=10;
      }
      for (int i=(mesh->ngz-1); i<2*(mesh->ngz-1); i++) {
	d_nnz[i]=20;
	d_nnz[localN-1-i]=20;
	o_nnz[i]=5;
	o_nnz[localN-1-i]=5;
      }
    }
    
    for (int i=2*(mesh->ngz-1); i<localN-2*((mesh->ngz-1));i++) {
      d_nnz[i]=25;
	d_nnz[localN-1-i]=25;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=0;
    }
    
    if (mesh->getNXPE()>1) {
      MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
    }
    else {
      MatSeqAIJSetPreallocation( MatA, 0, d_nnz );
    }
  }
  else {
    // first and last (mesh-ngz-1) entries are the edge x-values that (may) have 'off-diagonal' components (i.e. on another processor)
    if ( mesh->firstX() && mesh->lastX() ) {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=6;
	d_nnz[localN-1-i]=6;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=0;
      }
    }
    else if ( mesh->firstX() ) {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=6;
	d_nnz[localN-1-i]=6;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=3;
      }
    }
    else if ( mesh->lastX() ) {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=6;
	d_nnz[localN-1-i]=6;
	o_nnz[i]=3;
	o_nnz[localN-1-i]=0;
      }
    }
    else {
      for (int i=0; i<mesh->ngz-1; i++) {
	d_nnz[i]=6;
	d_nnz[localN-1-i]=6;
	o_nnz[i]=3;
	o_nnz[localN-1-i]=3;
      }
    }
    
    for (int i=mesh->ngz-1; i<localN-(mesh->ngz-1);i++) {
      d_nnz[i]=9;
	d_nnz[localN-1-i]=9;
	o_nnz[i]=0;
	o_nnz[localN-1-i]=0;
    }
    
    if (mesh->getNXPE()>1) {
      MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
    }
    else {
      MatSeqAIJSetPreallocation( MatA, 0, d_nnz );
    }
  }
  PetscFree( d_nnz );
  PetscFree( o_nnz );
  MatSetUp(MatA); 

  // Declare KSP Context 
  KSPCreate( comm, &ksp ); 

  // Get KSP Solver Type
  string type;
  opts->get("ksptype", type, KSP_GMRES);
  
  if(strcasecmp(type.c_str(), KSP_RICHARDSON) == 0)     ksptype = KSPRICHARDSON;
  else if(strcasecmp(type.c_str(), KSP_CHEBYSHEV) == 0) ksptype = KSPCHEBYSHEV;
  else if(strcasecmp(type.c_str(), KSP_CG) == 0)        ksptype = KSPCG;
  else if(type == "cgne")			ksptype = KSPCGNE;
  else if(type == "nash")			ksptype = KSPNASH;
  else if(type == "stcg")			ksptype = KSPSTCG;
  else if(type == "gltr")			ksptype = KSPGLTR;
  else if(strcasecmp(type.c_str(), KSP_GMRES) == 0)     ksptype = KSPGMRES;
  else if(type == "fgmres")			ksptype = KSPFGMRES;
  else if(type == "lgmres")			ksptype = KSPLGMRES;
  else if(type == "dgmres")			ksptype = KSPDGMRES;
  else if(type == "pgmres")			ksptype = KSPPGMRES;
  else if(strcasecmp(type.c_str(), KSP_TCQMR) == 0)     ksptype = KSPTCQMR;
  else if(strcasecmp(type.c_str(), KSP_BCGS) == 0)      ksptype = KSPBCGS;
  else if(type == "ibcgs")			ksptype = KSPIBCGS;
  else if(type == "fbcgs")			ksptype = KSPFBCGS;
  else if(type == "bcgsl")			ksptype = KSPBCGSL;
  else if(strcasecmp(type.c_str(), KSP_CGS) == 0)       ksptype = KSPCGS;
  else if(strcasecmp(type.c_str(), KSP_TFQMR) == 0)     ksptype = KSPTFQMR;
  else if(strcasecmp(type.c_str(), KSP_CR) == 0)        ksptype = KSPCR;
  else if(strcasecmp(type.c_str(), KSP_LSQR) == 0)      ksptype = KSPLSQR;
  else if(strcasecmp(type.c_str(), KSP_BICG) == 0)      ksptype = KSPBICG;
  else if(strcasecmp(type.c_str(), KSP_PREONLY) == 0)   ksptype = KSPPREONLY;
  else if(type == "qcg")			ksptype = KSPQCG;
  else if(type == "bicg")			ksptype = KSPBICG;
  else if(type == "minres")			ksptype = KSPMINRES;
  else if(type == "symmlq")			ksptype = KSPSYMMLQ;
  else if(type == "lcd")			ksptype = KSPLCD;
  else if(type == "python")			ksptype = KSPPYTHON;
  else if(type == "gcr")			ksptype = KSPGCR;
  else if(type == "specest")			ksptype = KSPSPECEST;
  else 
    throw BoutException("Unknown Krylov solver type '%s'", type.c_str());
  
  // Get preconditioner type
  // WARNING: only a few of these options actually make sense: see the PETSc documentation to work out which they are (possibly pbjacobi, sor might be useful choices?)
  string pctypeoption;
  opts->get("pctype", pctypeoption, "none", true);
  if (pctypeoption == "none") pctype = PCNONE;
  else if (pctypeoption == "jacobi") pctype = PCJACOBI;
  else if (pctypeoption == "sor") pctype = PCSOR;
  else if (pctypeoption == "lu") pctype = PCLU;
  else if (pctypeoption == "shell") pctype = PCSHELL;
  else if (pctypeoption == "bjacobi") pctype = PCBJACOBI;
  else if (pctypeoption == "mg") pctype = PCMG;
  else if (pctypeoption == "eisenstat") pctype = PCEISENSTAT;
  else if (pctypeoption == "ilu") pctype = PCILU;
  else if (pctypeoption == "icc") pctype = PCICC;
  else if (pctypeoption == "asm") pctype = PCASM;
  else if (pctypeoption == "gasm") pctype = PCGASM;
  else if (pctypeoption == "ksp") pctype = PCKSP;
  else if (pctypeoption == "composite") pctype = PCCOMPOSITE;
  else if (pctypeoption == "redundant") pctype = PCREDUNDANT;
  else if (pctypeoption == "spai") pctype = PCSPAI;
  else if (pctypeoption == "nn") pctype = PCNN;
  else if (pctypeoption == "cholesky") pctype = PCCHOLESKY;
  else if (pctypeoption == "pbjacobi") pctype = PCPBJACOBI;
  else if (pctypeoption == "mat") pctype = PCMAT;
  else if (pctypeoption == "hypre") pctype = PCHYPRE;
  else if (pctypeoption == "parms") pctype = PCPARMS;
  else if (pctypeoption == "fieldsplit") pctype = PCFIELDSPLIT;
  else if (pctypeoption == "tfs") pctype = PCTFS;
  else if (pctypeoption == "ml") pctype = PCML;
  else if (pctypeoption == "galerkin") pctype = PCGALERKIN;
  else if (pctypeoption == "exotic") pctype = PCEXOTIC;
  else if (pctypeoption == "hmpi") pctype = PCHMPI;
  else if (pctypeoption == "supportgraph") pctype = PCSUPPORTGRAPH;
  else if (pctypeoption == "asa") pctype = PCASA;
  else if (pctypeoption == "cp") pctype = PCCP;
  else if (pctypeoption == "bfbt") pctype = PCBFBT;
  else if (pctypeoption == "lsc") pctype = PCLSC;
  else if (pctypeoption == "python") pctype = PCPYTHON;
  else if (pctypeoption == "pfmg") pctype = PCPFMG;
  else if (pctypeoption == "syspfmg") pctype = PCSYSPFMG;
  else if (pctypeoption == "redistribute") pctype = PCREDISTRIBUTE;
  else if (pctypeoption == "svd") pctype = PCSVD;
  else if (pctypeoption == "gamg") pctype = PCGAMG;
  else if (pctypeoption == "sacusp") pctype = PCSACUSP;            /* these four run on NVIDIA GPUs using CUSP */
  else if (pctypeoption == "sacusppoly") pctype = PCSACUSPPOLY;
  else if (pctypeoption == "bicgstabcusp") pctype = PCBICGSTABCUSP;
  else if (pctypeoption == "ainvcusp") pctype = PCAINVCUSP;
  else if (pctypeoption == "bddc") pctype = PCBDDC;
  else 
    throw BoutException("Unknown KSP preconditioner type '%s'", pctypeoption.c_str());

  // Get Options specific to particular solver types
  opts->get("richardson_damping_factor",richardson_damping_factor,1.0,true);
  opts->get("chebyshev_max",chebyshev_max,100,true);
  opts->get("chebyshev_min",chebyshev_min,0.01,true);
  opts->get("gmres_max_steps",gmres_max_steps,30,true);
 
  // Get Tolerances for KSP solver
  opts->get("rtol",rtol,pow(10.0,-5),true);
  opts->get("atol",atol,pow(10.0,-50),true);
  opts->get("dtol",dtol,pow(10.0,5),true);
  opts->get("maxits",maxits,pow(10,5),true);

  // Get direct solver switch
  opts->get("direct", direct, false);
  if(direct)
    {
      output << endl << "Using LU decompostion for direct solution of system" << endl << endl;
    }
  
  // Ensure that the matrix is constructed first time
//   coefchanged = true;
  lastflag = -1;
}

const FieldPerp LaplacePetsc::solve(const FieldPerp &b) {
  return solve(b,b);
}

const FieldPerp LaplacePetsc::solve(const FieldPerp &b, const FieldPerp &x0) {
  #ifdef CHECK
    if ( flags & !implemented_flags) {
      if (flags&INVERT_4TH_ORDER) output<<"For PETSc based Laplacian inverter, use 'fourth_order=true' instead of setting INVERT_4TH_ORDER flag"<<endl;
      throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in petsc_laplace.cxx");
    }
  #endif
  
  int y = b.getIndex();           // Get the Y index
  sol.setIndex(y);// Initialize the solution field.
  sol = 0.;

  // Determine which row/columns of the matrix are locally owned
  MatGetOwnershipRange( MatA, &Istart, &Iend );

  int i = Istart;
{ Timer timer("petscsetup");
  
//     if ((fourth_order) && !(lastflag&INVERT_4TH_ORDER)) throw BoutException("Should not change INVERT_4TH_ORDER flag in LaplacePetsc: 2nd order and 4th order require different pre-allocation to optimize PETSc solver");
    
  // Set Matrix Elements

  // Loop over locally owned rows of matrix A - i labels NODE POINT from bottom left (0,0) = 0 to top right (meshx-1,meshz-1) = meshx*meshz-1
  // i increments by 1 for an increase of 1 in Z and by meshz for an increase of 1 in X.
  
  // X=0 to mesh->xstart-1 defines the boundary region of the domain.
  if( mesh->firstX() ) 
    {
      for(int x=0; x<mesh->xstart; x++)
	{
	  for(int z=0; z<mesh->ngz-1; z++) 
	    {
	      PetscScalar val;
	      // Set values corresponding to nodes adjacent in x if Neumann Boundary Conditions are required.
	      if(flags & INVERT_AC_IN_GRAD) 
		{
		  if( fourth_order )
		    {
		      // Fourth Order Accuracy on Boundary
		      Element(i,x,z, 0, 0, -25.0 / (12.0*mesh->dx[x][y]), MatA ); 
		      Element(i,x,z, 1, 0,   4.0 / mesh->dx[x][y], MatA ); 
		      Element(i,x,z, 2, 0,  -3.0 / mesh->dx[x][y], MatA );
		      Element(i,x,z, 3, 0,   4.0 / (3.0*mesh->dx[x][y]), MatA ); 
		      Element(i,x,z, 4, 0,  -1.0 / (4.0*mesh->dx[x][y]), MatA );
		    }
		  else
		    {
		      // Second Order Accuracy on Boundary
		      Element(i,x,z, 0, 0, -3.0 / (2.0*mesh->dx[x][y]), MatA ); 
		      Element(i,x,z, 1, 0,  2.0 / mesh->dx[x][y], MatA ); 
		      Element(i,x,z, 2, 0, -1.0 / (2.0*mesh->dx[x][y]), MatA ); 
// 			Element(i,x,z, 3, 0, 0.0, MatA );  // Reset these elements to 0 in case 4th order flag was used previously: not allowed now
// 			Element(i,x,z, 4, 0, 0.0, MatA );
		    }
		}
	      else
		{
		  // Set Diagonal Values to 1
		  val = 1;
		  Element(i,x,z, 0, 0, val, MatA );
		  
		  // Set off diagonal elements to zero
		  Element(i,x,z, 1, 0, 0.0, MatA );
		  Element(i,x,z, 2, 0, 0.0, MatA );
		  Element(i,x,z, 3, 0, 0.0, MatA );
		  Element(i,x,z, 4, 0, 0.0, MatA );
		}
	      
	      // Set Components of RHS and trial solution
	      val=0;
	      if( flags & INVERT_IN_RHS )         val = b[x][z];
	      else if( flags & INVERT_IN_SET )    val = x0[x][z];
	      VecSetValues( bs, 1, &i, &val, INSERT_VALUES );

	      val = x0[x][z];
	      VecSetValues( xs, 1, &i, &val, INSERT_VALUES );
	      
	      i++; // Increment row in Petsc matrix
	    }
	}
    }
  
  // Main domain with Laplacian operator
  for(int x=mesh->xstart; x <= mesh->xend; x++)
    {
      for(int z=0; z<mesh->ngz-1; z++) 
	{
	  BoutReal A0, A1, A2, A3, A4, A5;
	  A0 = A[x][y][z];
	  Coeffs( x, y, z, A1, A2, A3, A4, A5 );
	  
	  BoutReal dx   = mesh->dx[x][y];
	  BoutReal dx2  = pow( mesh->dx[x][y] , 2.0 );
	  BoutReal dz   = mesh->dz;
	  BoutReal dz2  = pow( mesh->dz, 2.0 );
	  BoutReal dxdz = mesh->dx[x][y] * mesh->dz;
	  
	  // Set Matrix Elements
	  PetscScalar val=0.;
	  if (fourth_order) {
	    // f(i,j) = f(x,z)
	    val = A0 - (5.0/2.0)*( (A1 / dx2) + (A2 / dz2) );
	    Element(i,x,z, 0, 0, val, MatA );
	    
	    // f(i-2,j-2)
	    val = A3 / ( 144.0 * dxdz );
	    Element(i,x,z, -2, -2, val, MatA ); 

	    // f(i-2,j-1)
	    val = -1.0 * A3 / ( 18.0 * dxdz );
	    Element(i,x,z, -2, -1, val, MatA );

	    // f(i-2,j)
	    val = (1.0/12.0) * ( (-1.0 * A1 /  dx2 ) + (A4 / dx) ); 
	    Element(i,x,z, -2, 0, val, MatA );

	    // f(i-2,j+1)
	    val = A3 / ( 18.0 * dxdz );
	    Element(i,x,z, -2, 1, val, MatA );

	    // f(i-2,j+2)
	    val = -1.0 * A3 / ( 144.0 * dxdz );
	    Element(i,x,z, -2, 2, val, MatA ); 

	    // f(i-1,j-2)
	    val = -1.0 * A3 / ( 18.0 * dxdz );
	    Element(i,x,z, -1, -2, val, MatA );

	    // f(i-1,j-1)
	    val = 4.0 * A3 / ( 9.0 * dxdz ); 
	    Element(i,x,z, -1, -1, val, MatA ); 

	    // f(i-1,j)
	    val = ( 4.0 * A1 / ( 3.0 * dx2 ) ) - ( 2.0 * A4 / ( 3.0 * dx ) );
	    Element(i,x,z, -1, 0, val, MatA );

	    // f(i-1,j+1)
	    val = -4.0 * A3 / ( 9.0 * dxdz ); 
	    Element(i,x,z, -1, 1, val, MatA );

	    // f(i-1,j+2)
	    val = A3 / ( 18.0 * dxdz );
	    Element(i,x,z, -1, 2, val, MatA );

	    // f(i,j-2)
	    val = (1.0/12.0) * ( ( -1.0 * A2 / dz2 ) + ( A5 / dz ) ); 
	    Element(i,x,z, 0, -2, val, MatA );
	    
	    // f(i,j-1)
	    val = ( 4.0 * A2 / ( 3.0 * dz2 ) ) - ( 2.0 * A5 / ( 3.0 * dz ) );
	    Element(i,x,z, 0, -1, val, MatA ); 

	    // f(i,j+1)
	    val = ( 4.0 * A2 / ( 3.0 * dz2 ) ) + ( 2.0 * A5 / ( 3.0 * dz ) );
	    Element(i,x,z, 0, 1, val, MatA );

	    // f(i,j+2)
	    val = (-1.0/12.0) * ( ( A2 / dz2 ) + ( A5 / dz ) ); 
	    Element(i,x,z, 0, 2, val, MatA );

	    // f(i+1,j-2)
	    val = A3 / ( 18.0 * dxdz );
	    Element(i,x,z, 1, -2, val, MatA );
	    
	    // f(i+1,j-1)
	    val = -4.0 * A3 / ( 9.0 * dxdz );
	    Element(i,x,z, 1, -1, val, MatA );
	    
	    // f(i+1,j)
	    val = ( 4.0 * A1 / ( 3.0*dx2 ) ) + ( 2.0 * A4 / ( 3.0 * dx ) ); 
	    Element(i,x,z, 1, 0, val, MatA );
	    
	    // f(i+1,j+1)
	    val = 4.0 * A3 / ( 9.0 * dxdz ); 
	    Element(i,x,z, 1, 1, val, MatA );

	    // f(i+1,j+2)
	    val = -1.0 * A3 / ( 18.0 * dxdz );
	    Element(i,x,z, 1, 2, val, MatA );

	    // f(i+2,j-2)
	    val = -1.0 * A3 / ( 144.0 * dxdz );
	    Element(i,x,z, 2, -2, val, MatA ); 

	    // f(i+2,j-1)
	    val = A3 / ( 18.0 * dxdz );
	    Element(i,x,z, 2, -1, val, MatA );

	    // f(i+2,j)
	    val = (-1.0/12.0) * ( (A1 / dx2) + (A4 / dx) ); 
	    Element(i,x,z, 2, 0, val, MatA );

	    // f(i+2,j+1)
	    val = -1.0 * A3 / ( 18.0 * dxdz );
	    Element(i,x,z, 2, 1, val, MatA );

	    // f(i+2,j+2)
	    val = A3 / ( 144.0 * dxdz );
	    Element(i,x,z, 2, 2, val, MatA ); 
	  }
	  else {
	    // f(i,j) = f(x,z)
	    val = A0 - 2.0*( (A1 / dx2) + (A2 / dz2) ); 
	    Element(i,x,z, 0, 0, val, MatA );
	    
	    // f(i-1,j-1)
	    val = A3 / (4.0 * dxdz); 
	    Element(i,x,z, -1, -1, val, MatA ); 

	    // f(i-1,j)
	    val = ( A1 / dx2 ) - A4 / ( 2.0 * dx ); 
	    Element(i,x,z, -1, 0, val, MatA );

	    // f(i-1,j+1)
	    val = -1.0 * A3 / ( 4.0 * dxdz ); 
	    Element(i,x,z, -1, 1, val, MatA );

	    // f(i,j-1)
	    val = ( A2 / dz2 ) - ( A5 / ( 2.0 * dz ) ); 
	    Element(i,x,z, 0, -1, val, MatA ); 

	    // f(i,j+1)
	    val = ( A2 / dz2 ) + ( A5 / ( 2.0 * dz ) ); 
	    Element(i,x,z, 0, 1, val, MatA );

	    // f(i+1,j-1)
	    val = -1.0 * A3 / ( 4.0 * dxdz ); 
	    Element(i,x,z, 1, -1, val, MatA );
	    
	    // f(i+1,j)
	    val = ( A1 / dx2 ) + ( A4 / ( 2.0 * dx ) ); 
	    Element(i,x,z, 1, 0, val, MatA );
	    
	    // f(i+1,j+1)
	    val = A3 / ( 4.0 * dxdz ); 
	    Element(i,x,z, 1, 1, val, MatA );
	  }
	  // Set Components of RHS Vector
	  val  = b[x][z];
	  VecSetValues( bs, 1, &i, &val, INSERT_VALUES );

	  // Set Components of Trial Solution Vector
	  val = x0[x][z];
	  VecSetValues( xs, 1, &i, &val, INSERT_VALUES ); 
	  i++;
	}
    }
  
  // X=mesh->xend+1 to mesh->ngx-1 defines the upper boundary region of the domain.
  if( mesh->lastX() ) 
    {
      for(int x=mesh->xend+1; x<mesh->ngx; x++)
	{
	  for(int z=0; z<mesh->ngz-1; z++) 
	    {
	      // Set Diagonal Values to 1
	      PetscScalar val = 1;
	      Element(i,x,z, 0, 0, val, MatA ); 
	      
	      // Set values corresponding to nodes adjacent in x if Neumann Boundary Conditions are required.
	      if(flags & INVERT_AC_OUT_GRAD) 
		{
		  if( fourth_order )
		    {
		      // Fourth Order Accuracy on Boundary
		      Element(i,x,z,  0, 0, 25.0 / (12.0*mesh->dx[x][y]), MatA ); 
		      Element(i,x,z, -1, 0, -4.0 / mesh->dx[x][y], MatA ); 
		      Element(i,x,z, -2, 0,  3.0 / mesh->dx[x][y], MatA );
		      Element(i,x,z, -3, 0, -4.0 / (3.0*mesh->dx[x][y]), MatA ); 
		      Element(i,x,z, -4, 0,  1.0 / (4.0*mesh->dx[x][y]), MatA );
		    }
		  else
		    {
		      // Second Order Accuracy on Boundary
		      Element(i,x,z,  0, 0,  3.0 / (2.0*mesh->dx[x][y]), MatA ); 
		      Element(i,x,z, -1, 0, -2.0 / mesh->dx[x][y], MatA ); 
		      Element(i,x,z, -2, 0,  1.0 / (2.0*mesh->dx[x][y]), MatA ); 
// 			Element(i,x,z, -3, 0,  0.0, MatA );  // Reset these elements to 0 in case 4th order flag was used previously: not allowed now
// 			Element(i,x,z, -4, 0,  0.0, MatA );
		    }
		}
	      else
		{
		  // Set off diagonal elements to zero
		  Element(i,x,z, -1, 0, 0.0, MatA );
		  Element(i,x,z, -2, 0, 0.0, MatA );
		  Element(i,x,z, -3, 0, 0.0, MatA );
		  Element(i,x,z, -4, 0, 0.0, MatA );
		}            
	      
	      // Set Components of RHS
	      val=0;
	      if( flags & INVERT_OUT_RHS )        val = b[x][z];
	      else if( flags & INVERT_OUT_SET )   val = x0[x][z];
	      VecSetValues( bs, 1, &i, &val, INSERT_VALUES );

	      val = x0[x][z];
	      VecSetValues( xs, 1, &i, &val, INSERT_VALUES );
	      
	      i++; // Increment row in Petsc matrix
	    }
	}
    }
  
  if(i != Iend) {
    throw BoutException("Petsc index sanity check failed");
  }
  
  // Assemble Matrix
  MatAssemblyBegin( MatA, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatA, MAT_FINAL_ASSEMBLY );

  // Record which flags were used for this matrix
  lastflag = flags;
  
  // Assemble RHS Vector
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);

  // Assemble Trial Solution Vector
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);

  // Configure Linear Solver               
  KSPSetOperators( ksp,MatA,MatA,DIFFERENT_NONZERO_PATTERN );

  PC pc;
 
  if(direct)
    {
      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      PCFactorSetMatSolverPackage(pc,"mumps");
    }
  else
    {
      KSPSetType( ksp, ksptype );
      
      if( ksptype == KSPRICHARDSON )     KSPRichardsonSetScale( ksp, richardson_damping_factor );
      else if( ksptype == KSPCHEBYSHEV ) KSPChebyshevSetEigenvalues( ksp, chebyshev_max, chebyshev_min );
      else if( ksptype == KSPGMRES )     KSPGMRESSetRestart( ksp, gmres_max_steps );
      
      KSPSetTolerances( ksp, rtol, atol, dtol, maxits );
      
      if( !( flags & INVERT_START_NEW ) ) KSPSetInitialGuessNonzero( ksp, (PetscBool) true );
      
      KSPGetPC(ksp,&pc);
      PCSetType(pc, pctype);
      
      KSPSetFromOptions( ksp );           
    }
  
}
{ Timer timer("petscsolve");
  // Solve the system
  KSPSolve( ksp, bs, xs );
}
  
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );
  if (reason==-3) { // Too many iterations, might be fixed by taking smaller timestep
    throw BoutIterationFail("petsc_laplace: too many iterations");
  }
  else if (reason<=0) {
    output<<"KSPConvergedReason is "<<reason<<endl;
    throw BoutException("petsc_laplace: inversion failed to converge.");
  }

  // Add data to FieldPerp Object
  i = Istart;
  if(mesh->firstX()) 
    {
      for(int x=0; x<mesh->xstart; x++)
	{
	  for(int z=0; z<mesh->ngz-1; z++) 
	    {
	      PetscScalar val = 0;
	      VecGetValues(xs, 1, &i, &val ); 
	      sol[x][z] = val;
	      i++; // Increment row in Petsc matrix
	    }
	}
    }
  
  for(int x=mesh->xstart; x <= mesh->xend; x++)
    {
      for(int z=0; z<mesh->ngz-1; z++) 
	{
	  PetscScalar val = 0;
	  VecGetValues(xs, 1, &i, &val ); 
	  sol[x][z] = val;
	  i++; // Increment row in Petsc matrix
	}
    }

  if(mesh->lastX()) 
    {
      for(int x=mesh->xend+1; x<mesh->ngx; x++)
	{
	  for(int z=0;z < mesh->ngz-1; z++) 
	    {
	      PetscScalar val = 0;
	      VecGetValues(xs, 1, &i, &val ); 
	      sol[x][z] = val;
	      i++; // Increment row in Petsc matrix
	    }	
	}
    }
  
  if(i != Iend) {
    throw BoutException("Petsc index sanity check 2 failed");
  }

  return sol;
}

void LaplacePetsc::Element(int i, int x, int z, int xshift, int zshift, PetscScalar ele, Mat &MatA ) 
{
  // Need to convert LOCAL x to GLOBAL x in order to correctly calculate PETSC Matrix Index. D'oh!

  int xoffset = Istart / meshz;
  if( Istart % meshz != 0 )
    throw  BoutException("Petsc index sanity check 3 failed");

  int row_new = x + xshift; // should never be out of range.
  if( !mesh->firstX() ) row_new += (xoffset - mesh->xstart);

  int col_new = z + zshift;
  if( col_new < 0 )            col_new += meshz;
  else if( col_new > meshz-1 ) col_new -= meshz;
 
  // convert to global indices
  int index = (row_new * meshz) + col_new;
  
  MatSetValues(MatA,1,&i,1,&index,&ele,INSERT_VALUES); 
}

void LaplacePetsc::Coeffs( int x, int y, int z, BoutReal &coef1, BoutReal &coef2, BoutReal &coef3, BoutReal &coef4, BoutReal &coef5 )
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
//   if( (x > 0) && (x < (mesh->ngx-1)) ) //Valid if doing second order derivative, not if fourth: should only be called for xstart<=x<=xend anyway
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

#endif // BOUT_HAS_PETSC_3_3
