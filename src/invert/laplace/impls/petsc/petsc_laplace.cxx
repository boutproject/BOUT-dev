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
#ifdef BOUT_HAS_PETSC

#include "petsc_laplace.hxx"

#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <bout/assert.hxx>
#include <utils.hxx>

#define KSP_RICHARDSON "richardson"
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

#undef __FUNCT__
#define __FUNCT__ "laplacePCapply"
static PetscErrorCode laplacePCapply(PC pc,Vec x,Vec y) {
  int ierr;

  // Get the context
  LaplacePetsc *s;
  ierr = PCShellGetContext(pc,(void**)&s);CHKERRQ(ierr);

  PetscFunctionReturn(s->precon(x, y));
}

LaplacePetsc::LaplacePetsc(Options *opt, const CELL_LOC loc, Mesh *mesh_in) :
  Laplacian(opt, loc, mesh_in),
  A(0.0), C1(1.0), C2(1.0), D(1.0), Ex(0.0), Ez(0.0),
  issetD(false), issetC(false), issetE(false),
  lib(opt==nullptr ? &(Options::root()["laplace"]["petsc"]) : &(*opt)["petsc"])
{
  A.setLocation(location);
  C1.setLocation(location);
  C2.setLocation(location);
  D.setLocation(location);
  Ex.setLocation(location);
  Ez.setLocation(location);

  // Get Options in Laplace Section
  if (!opt) opts = Options::getRoot()->getSection("laplace");
  else opts=opt;

  #if CHECK > 0
    // These are the implemented flags
    implemented_flags = INVERT_START_NEW;
    implemented_boundary_flags = INVERT_AC_GRAD
                                 + INVERT_SET
                                 + INVERT_RHS
                                 ;
    // Checking flags are set to something which is not implemented
    // This is done binary (which is possible as each flag is a power of 2)
    if ( global_flags & ~implemented_flags ) {
      if (global_flags&INVERT_4TH_ORDER) output<<"For PETSc based Laplacian inverter, use 'fourth_order=true' instead of setting INVERT_4TH_ORDER flag"<<endl;
      throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in petsc_laplace.cxx");
    }
    if ( inner_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if ( outer_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if(localmesh->periodicX) {
      throw BoutException("LaplacePetsc does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
      }
  #endif

  // Get communicator for group of processors in X - all points in z-x plane for fixed y.
  comm = localmesh->getXcomm();

  // Need to determine local size to use based on prior parallelisation
  // Coefficient values are stored only on local processors.
  localN = (localmesh->xend - localmesh->xstart + 1) * (localmesh->LocalNz);
  if(localmesh->firstX())
    localN += localmesh->xstart * (localmesh->LocalNz);    // If on first processor add on width of boundary region
  if(localmesh->lastX())
    localN += localmesh->xstart * (localmesh->LocalNz);    // If on last processor add on width of boundary region

  // Calculate 'size' (the total number of points in physical grid)
  if(MPI_Allreduce(&localN, &size, 1, MPI_INT, MPI_SUM, comm) != MPI_SUCCESS)
    throw BoutException("Error in MPI_Allreduce during LaplacePetsc initialisation");

  // Calculate total (physical) grid dimensions
  meshz = localmesh->LocalNz;
  meshx = size / meshz;

  // Create PETSc type of vectors for the solution and the RHS vector
  VecCreate( comm, &xs );
  VecSetSizes( xs, localN, size );
  VecSetFromOptions( xs );
  VecDuplicate( xs , &bs );

  // Get 4th order solver switch
  opts->get("fourth_order", fourth_order, false);

  // Set size of (the PETSc) Matrix on each processor to localN x localN
  MatCreate( comm, &MatA );
  MatSetSizes( MatA, localN, localN, size, size );
  MatSetFromOptions(MatA);
  //   if (fourth_order) MatMPIAIJSetPreallocation( MatA, 25, PETSC_NULL, 10, PETSC_NULL );
//   else MatMPIAIJSetPreallocation( MatA, 9, PETSC_NULL, 3, PETSC_NULL );

  /* Pre allocate memory
   * nnz denotes an array containing the number of non-zeros in the various rows
   * for
   * d_nnz - The diagonal terms in the matrix
   * o_nnz - The off-diagonal terms in the matrix (needed when running in
   *         parallel)
   */
  PetscInt *d_nnz, *o_nnz;
  PetscMalloc( (localN)*sizeof(PetscInt), &d_nnz );
  PetscMalloc( (localN)*sizeof(PetscInt), &o_nnz );
  if (fourth_order) {
    // first and last 2*localmesh-LocalNz entries are the edge x-values that (may) have 'off-diagonal' components (i.e. on another processor)
    if ( localmesh->firstX() && localmesh->lastX() ) {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=15;
        d_nnz[localN-1-i]=15;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=0;
      }
      for (int i=(localmesh->LocalNz); i<2*(localmesh->LocalNz); i++) {
        d_nnz[i]=20;
        d_nnz[localN-1-i]=20;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=0;
      }
    }
    else if ( localmesh->firstX() ) {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=15;
        d_nnz[localN-1-i]=15;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=10;
      }
      for (int i=(localmesh->LocalNz); i<2*(localmesh->LocalNz); i++) {
        d_nnz[i]=20;
        d_nnz[localN-1-i]=20;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=5;
      }
    }
    else if ( localmesh->lastX() ) {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=15;
        d_nnz[localN-1-i]=15;
        o_nnz[i]=10;
        o_nnz[localN-1-i]=0;
      }
      for (int i=(localmesh->LocalNz); i<2*(localmesh->LocalNz); i++) {
        d_nnz[i]=20;
        d_nnz[localN-1-i]=20;
        o_nnz[i]=5;
        o_nnz[localN-1-i]=0;
      }
    }
    else {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=15;
        d_nnz[localN-1-i]=15;
        o_nnz[i]=10;
        o_nnz[localN-1-i]=10;
      }
      for (int i=(localmesh->LocalNz); i<2*(localmesh->LocalNz); i++) {
        d_nnz[i]=20;
        d_nnz[localN-1-i]=20;
        o_nnz[i]=5;
        o_nnz[localN-1-i]=5;
      }
    }

    for (int i=2*(localmesh->LocalNz); i<localN-2*((localmesh->LocalNz));i++) {
      d_nnz[i]=25;
        d_nnz[localN-1-i]=25;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=0;
    }

    // Use d_nnz and o_nnz for preallocating the matrix
    if (localmesh->firstX() && localmesh->lastX()) {
      // Only one processor in X
      MatSeqAIJSetPreallocation( MatA, 0, d_nnz );
    }else {
      MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
    }
  }
  else {
    // first and last localmesh->LocalNz entries are the edge x-values that (may) have 'off-diagonal' components (i.e. on another processor)
    if ( localmesh->firstX() && localmesh->lastX() ) {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=6;
        d_nnz[localN-1-i]=6;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=0;
      }
    }
    else if ( localmesh->firstX() ) {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=6;
        d_nnz[localN-1-i]=6;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=3;
      }
    }
    else if ( localmesh->lastX() ) {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=6;
        d_nnz[localN-1-i]=6;
        o_nnz[i]=3;
        o_nnz[localN-1-i]=0;
      }
    }
    else {
      for (int i=0; i<localmesh->LocalNz; i++) {
        d_nnz[i]=6;
        d_nnz[localN-1-i]=6;
        o_nnz[i]=3;
        o_nnz[localN-1-i]=3;
      }
    }

    for (int i=localmesh->LocalNz; i<localN-(localmesh->LocalNz);i++) {
      d_nnz[i]=9;
        d_nnz[localN-1-i]=9;
        o_nnz[i]=0;
        o_nnz[localN-1-i]=0;
    }

    // Use d_nnz and o_nnz for preallocating the matrix
    if (localmesh->firstX() && localmesh->lastX()) {
      MatSeqAIJSetPreallocation( MatA, 0, d_nnz );
    } else {
      MatMPIAIJSetPreallocation( MatA, 0, d_nnz, 0, o_nnz );
    }
  }
  // Free the d_nnz and o_nnz arrays, as these are will not be used anymore
  PetscFree( d_nnz );
  PetscFree( o_nnz );
  // Sets up the internal matrix data structures for the later use.
  MatSetUp(MatA);

  // Declare KSP Context (abstract PETSc object that manages all Krylov methods)
  ksp = lib.createKSPWithOptions(comm);

  // Get KSP Solver Type (Generalizes Minimal RESidual is the default)
  ksptype = (*opts)["ksptype"].doc("KSP solver type").withDefault(KSP_GMRES);

  // Get preconditioner type
  // WARNING: only a few of these options actually make sense: see the
  // PETSc documentation to work out which they are (possibly
  // pbjacobi, sor might be useful choices?)
  pctype = (*opts)["pctype"]
               .doc("Preconditioner type. See the PETSc documentation for options")
               .withDefault("none");

  // Let "user" be a synonym for "shell"
  if (pctype == "user") {
    pctype = PCSHELL;
  }
  
  // Get Options specific to particular solver types
  opts->get("richardson_damping_factor",richardson_damping_factor,1.0,true);
  opts->get("chebyshev_max",chebyshev_max,100,true);
  opts->get("chebyshev_min",chebyshev_min,0.01,true);
  opts->get("gmres_max_steps",gmres_max_steps,30,true);

  // Get Tolerances for KSP solver
  rtol = (*opts)["rtol"].doc("Relative tolerance for KSP solver").withDefault(1e-5);
  atol = (*opts)["atol"].doc("Absolute tolerance for KSP solver").withDefault(1e-50);
  dtol = (*opts)["dtol"].doc("Divergence tolerance for KSP solver").withDefault(1e5);
  maxits = (*opts)["maxits"].doc("Maximum number of KSP iterations").withDefault(100000);

  // Get direct solver switch
  direct = (*opts)["direct"].doc("Use direct (LU) solver?").withDefault(false);
  if (direct) {
    output << endl << "Using LU decompostion for direct solution of system" << endl << endl;
  }

  pcsolve = nullptr;
  if (pctype == PCSHELL) {

    rightprec = (*opts)["rightprec"].doc("Right preconditioning?").withDefault(true);

    // Options for preconditioner are in a subsection
    pcsolve = Laplacian::create(opts->getSection("precon"));
  }

  // Ensure that the matrix is constructed first time
  //   coefchanged = true;
  //  lastflag = -1;
}

FieldPerp LaplacePetsc::solve(const FieldPerp& b) { return solve(b, b); }

/*!
 * Solves Ax=b for x given a b and an initial guess for x (x0)
 *
 * This function will:
 *      1. Set the matrix element of the matrix A, used to solve Ax=b
 *         (this includes setting the values for the bounary condition)
 *      2. Solve the matrix Ax = b
 *
 * \param[in] b     The RHS of the equation Ax=b.
 *                  This is an y-slice of the original field. The field wil be
 *                  flattened to an 1D array in order to write the equation on
 *                  the form Ax=b
 * \param[in] x0    The initial guess for the solver.
 *                  May also contain the  boundary condition if flag 32 - INVERT_SET is set
 *
 * \returns sol     The solution x of the problem Ax=b.
 */
FieldPerp LaplacePetsc::solve(const FieldPerp& b, const FieldPerp& x0) {
  TRACE("LaplacePetsc::solve");

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  
  #if CHECK > 0
    // Checking flags are set to something which is not implemented (see
    // constructor for details)
    if ( global_flags & !implemented_flags) {
      if (global_flags&INVERT_4TH_ORDER) output<<"For PETSc based Laplacian inverter, use 'fourth_order=true' instead of setting INVERT_4TH_ORDER flag"<<endl;
      throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in petsc_laplace.cxx");
    }
    if ( inner_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
    if ( outer_boundary_flags & ~implemented_boundary_flags ) {
      throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in petsc_laplace.cxx");
    }
  #endif

  int y = b.getIndex(); // Get the Y index
  sol.setIndex(y);      // Initialize the solution field.
  sol = 0.;

  // Determine which row/columns of the matrix are locally owned
  MatGetOwnershipRange( MatA, &Istart, &Iend );

  int i = Istart;   // The row in the PETSc matrix
  { Timer timer("petscsetup");

    //     if ((fourth_order) && !(lastflag&INVERT_4TH_ORDER)) throw BoutException("Should not change INVERT_4TH_ORDER flag in LaplacePetsc: 2nd order and 4th order require different pre-allocation to optimize PETSc solver");

  /* Set Matrix Elements
   *
   * Loop over locally owned rows of matrix A
   * i labels NODE POINT from
   * bottom left = (0,0) = 0
   * to
   * top right = (meshx-1,meshz-1) = meshx*meshz-1
   *
   * i increments by 1 for an increase of 1 in Z
   * i increments by meshz for an increase of 1 in X.
   *
   * In other word the indexing is done in a row-major order, but starting at
   * bottom left rather than top left
   */
  // X=0 to localmesh->xstart-1 defines the boundary region of the domain.
  // Set the values for the inner boundary region
  if( localmesh->firstX() ) {
      for(int x=0; x<localmesh->xstart; x++) {
          for(int z=0; z<localmesh->LocalNz; z++) {
              PetscScalar val; // Value of element to be set in the matrix
              // If Neumann Boundary Conditions are set.
              if(inner_boundary_flags & INVERT_AC_GRAD) {
                  // Set values corresponding to nodes adjacent in x
                  if( fourth_order ) {
                      // Fourth Order Accuracy on Boundary
                      Element(i,x,z, 0, 0, -25.0 / (12.0*coords->dx(x,y)) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, 1, 0,   4.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, 2, 0,  -3.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, 3, 0,   4.0 / (3.0*coords->dx(x,y)) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, 4, 0,  -1.0 / (4.0*coords->dx(x,y)) / sqrt(coords->g_11(x,y)), MatA );
                  } else {
//                    // Second Order Accuracy on Boundary
//                    Element(i,x,z, 0, 0, -3.0 / (2.0*coords->dx(x,y)), MatA );
//                    Element(i,x,z, 1, 0,  2.0 / coords->dx(x,y), MatA );
//                    Element(i,x,z, 2, 0, -1.0 / (2.0*coords->dx(x,y)), MatA );
// //                   Element(i,x,z, 3, 0, 0.0, MatA );  // Reset these elements to 0 in case 4th order flag was used previously: not allowed now
// //                   Element(i,x,z, 4, 0, 0.0, MatA );
                      // Second Order Accuracy on Boundary, set half-way between grid points
                      Element(i,x,z, 0, 0, -1.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, 1, 0,  1.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, 2, 0, 0.0, MatA );
//                      Element(i,x,z, 3, 0, 0.0, MatA );  // Reset these elements to 0 in case 4th order flag was used previously: not allowed now
//                      Element(i,x,z, 4, 0, 0.0, MatA );
                    }
                } else {
                  if (fourth_order) {
                      // Set Diagonal Values to 1
                      Element(i,x,z, 0, 0, 1., MatA );

                      // Set off diagonal elements to zero
                      Element(i,x,z, 1, 0, 0.0, MatA );
                      Element(i,x,z, 2, 0, 0.0, MatA );
                      Element(i,x,z, 3, 0, 0.0, MatA );
                      Element(i,x,z, 4, 0, 0.0, MatA );
                    } else {
                      Element(i,x,z, 0, 0, 0.5, MatA );
                      Element(i,x,z, 1, 0, 0.5, MatA );
                      Element(i,x,z, 2, 0, 0., MatA );
                    }
                }

              val=0; // Initialize val

              // Set Components of RHS
              // If the inner boundary value should be set by b or x0
              if( inner_boundary_flags & INVERT_RHS )         val = b[x][z];
              else if( inner_boundary_flags & INVERT_SET )    val = x0[x][z];

              // Set components of the RHS (the PETSc vector bs)
              // 1 element is being set in row i to val
              // INSERT_VALUES replaces existing entries with new values
              VecSetValues( bs, 1, &i, &val, INSERT_VALUES );

              // Set components of the and trial solution (the PETSc vector xs)
              // 1 element is being set in row i to val
              // INSERT_VALUES replaces existing entries with new values
              val = x0[x][z];
              VecSetValues( xs, 1, &i, &val, INSERT_VALUES );

              i++; // Increment row in Petsc matrix
            }
        }
    }

  // Set the values for the main domain
  for(int x=localmesh->xstart; x <= localmesh->xend; x++) {
    for(int z=0; z<localmesh->LocalNz; z++) {
        // NOTE: Only A0 is the A from setCoefA ()
        BoutReal A0, A1, A2, A3, A4, A5;
        A0 = A(x,y,z);

        ASSERT3(finite(A0));
        
        // Set the matrix coefficients
        Coeffs( x, y, z, A1, A2, A3, A4, A5 );

        BoutReal dx   = coords->dx(x,y);
        BoutReal dx2  = SQ(coords->dx(x,y));
        BoutReal dz   = coords->dz;
        BoutReal dz2  = SQ(coords->dz);
        BoutReal dxdz = coords->dx(x,y) * coords->dz;
        
        ASSERT3(finite(A1));
        ASSERT3(finite(A2));
        ASSERT3(finite(A3));
        ASSERT3(finite(A4));
        ASSERT3(finite(A5));
        
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
          } else {
            // Second order
            
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

  // X=localmesh->xend+1 to localmesh->LocalNx-1 defines the upper boundary region of the domain.
  // Set the values for the outer boundary region
  if( localmesh->lastX() ) {
      for(int x=localmesh->xend+1; x<localmesh->LocalNx; x++) {
          for(int z=0; z<localmesh->LocalNz; z++) {
              // Set Diagonal Values to 1
              PetscScalar val = 1;
              Element(i,x,z, 0, 0, val, MatA );

              // If Neumann Boundary Conditions are set.
              if(outer_boundary_flags & INVERT_AC_GRAD) {
                  // Set values corresponding to nodes adjacent in x
                  if( fourth_order ) {
                      // Fourth Order Accuracy on Boundary
                      Element(i,x,z,  0, 0, 25.0 / (12.0*coords->dx(x,y)) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, -1, 0, -4.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, -2, 0,  3.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, -3, 0, -4.0 / (3.0*coords->dx(x,y)) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, -4, 0,  1.0 / (4.0*coords->dx(x,y)) / sqrt(coords->g_11(x,y)), MatA );
                    }
                  else {
//                    // Second Order Accuracy on Boundary
//                    Element(i,x,z,  0, 0,  3.0 / (2.0*coords->dx(x,y)), MatA );
//                    Element(i,x,z, -1, 0, -2.0 / coords->dx(x,y), MatA );
//                    Element(i,x,z, -2, 0,  1.0 / (2.0*coords->dx(x,y)), MatA );
// //                   Element(i,x,z, -3, 0,  0.0, MatA );  // Reset these elements to 0 in case 4th order flag was used previously: not allowed now
// //                   Element(i,x,z, -4, 0,  0.0, MatA );
                      // Second Order Accuracy on Boundary, set half-way between grid points
                      Element(i,x,z,  0, 0,  1.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, -1, 0, -1.0 / coords->dx(x,y) / sqrt(coords->g_11(x,y)), MatA );
                      Element(i,x,z, -2, 0,  0.0, MatA );
//                      Element(i,x,z, -3, 0,  0.0, MatA );  // Reset these elements to 0 in case 4th order flag was used previously: not allowed now
//                      Element(i,x,z, -4, 0,  0.0, MatA );
                    }
                }
              else {
                if (fourth_order) {
                    // Set off diagonal elements to zero
                    Element(i,x,z, -1, 0, 0.0, MatA );
                    Element(i,x,z, -2, 0, 0.0, MatA );
                    Element(i,x,z, -3, 0, 0.0, MatA );
                    Element(i,x,z, -4, 0, 0.0, MatA );
                  }
                else {
                    Element(i,x,z,  0, 0, 0.5 , MatA );
                    Element(i,x,z, -1, 0, 0.5 , MatA );
                    Element(i,x,z, -2, 0, 0., MatA );
                  }
              }

              // Set Components of RHS
              // If the inner boundary value should be set by b or x0
              val=0;
              if( outer_boundary_flags & INVERT_RHS )        val = b[x][z];
              else if( outer_boundary_flags & INVERT_SET )   val = x0[x][z];

              // Set components of the RHS (the PETSc vector bs)
              // 1 element is being set in row i to val
              // INSERT_VALUES replaces existing entries with new values
              VecSetValues( bs, 1, &i, &val, INSERT_VALUES );

              // Set components of the and trial solution (the PETSc vector xs)
              // 1 element is being set in row i to val
              // INSERT_VALUES replaces existing entries with new values
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

//   // Record which flags were used for this matrix
//   lastflag = flags;

  // Assemble RHS Vector
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);

  // Assemble Trial Solution Vector
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);

  // Configure Linear Solver
#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators( ksp,MatA,MatA);
#else
  KSPSetOperators( ksp,MatA,MatA,DIFFERENT_NONZERO_PATTERN );
#endif
  PC pc; // The preconditioner option

  if(direct) { // If a direct solver has been chosen
    // Get the preconditioner
    KSPGetPC(ksp,&pc);
    // Set the preconditioner
    PCSetType(pc,PCLU);
    // Set the solver type
#if PETSC_VERSION_GE(3,9,0)
    PCFactorSetMatSolverType(pc,"mumps");
#else
    PCFactorSetMatSolverPackage(pc,"mumps");
#endif
  }else { // If a iterative solver has been chosen
    KSPSetType( ksp, ksptype.c_str() ); // Set the type of the solver

    if( ksptype == KSPRICHARDSON )     KSPRichardsonSetScale( ksp, richardson_damping_factor );
#ifdef KSPCHEBYSHEV
    else if( ksptype == KSPCHEBYSHEV ) KSPChebyshevSetEigenvalues( ksp, chebyshev_max, chebyshev_min );
#endif
    else if( ksptype == KSPGMRES )     KSPGMRESSetRestart( ksp, gmres_max_steps );

    // Set the relative and absolute tolerances
    KSPSetTolerances( ksp, rtol, atol, dtol, maxits );

    // If the initial guess is not set to zero
    if( !( global_flags & INVERT_START_NEW ) ) KSPSetInitialGuessNonzero( ksp, (PetscBool) true );

    // Get the preconditioner
    KSPGetPC(ksp,&pc);

    // Set the type of the preconditioner
    PCSetType(pc, pctype.c_str());

    // If pctype = user in BOUT.inp, it will be translated to PCSHELL upon
    // construction of the object
    if (pctype == PCSHELL) {
      // User-supplied preconditioner function
      PCShellSetApply(pc,laplacePCapply);
      PCShellSetContext(pc,this);
      if(rightprec) {
        KSPSetPCSide(ksp, PC_RIGHT); // Right preconditioning
      }else
        KSPSetPCSide(ksp, PC_LEFT);  // Left preconditioning
      //ierr = PCShellSetApply(pc,laplacePCapply);CHKERRQ(ierr);
      //ierr = PCShellSetContext(pc,this);CHKERRQ(ierr);
      //ierr = KSPSetPCSide(ksp, PC_RIGHT);CHKERRQ(ierr);
    }

    KSPSetFromOptions( ksp );
  }
  }

  // Call the actual solver
  { Timer timer("petscsolve");
    KSPSolve( ksp, bs, xs ); // Call the solver to solve the system
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
  // Set the inner boundary values
  if(localmesh->firstX()) {
      for(int x=0; x<localmesh->xstart; x++) {
          for(int z=0; z<localmesh->LocalNz; z++) {
              PetscScalar val = 0;
              VecGetValues(xs, 1, &i, &val );
              sol[x][z] = val;
              i++; // Increment row in Petsc matrix
            }
        }
    }

  // Set the main domain values
  for(int x=localmesh->xstart; x <= localmesh->xend; x++) {
      for(int z=0; z<localmesh->LocalNz; z++) {
          PetscScalar val = 0;
          VecGetValues(xs, 1, &i, &val );
          sol[x][z] = val;
          i++; // Increment row in Petsc matrix
        }
    }

  // Set the outer boundary values
  if(localmesh->lastX()) {
      for(int x=localmesh->xend+1; x<localmesh->LocalNx; x++) {
          for(int z=0;z < localmesh->LocalNz; z++) {
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

  // Return the solution
  return sol;
}

/*!
 * Sets the elements of the matrix A, which is used to solve the problem Ax=b.
 *
 * \param[in]
 * i
 * The row of the PETSc matrix
 * \param[in] x         Local x index of the mesh
 * \param[in] z         Local z index of the mesh
 * \param[in] xshift    The shift in rows from the index x
 * \param[in] zshift    The shift in columns from the index z
 * \param[in] ele       Value of the element
 * \param[in] MatA      The matrix A used in the inversion
 *
 * \param[out] MatA     The matrix A used in the inversion
 */
void LaplacePetsc::Element(int i, int x, int z,
                           int xshift, int zshift,
                           PetscScalar ele, Mat &MatA ) {

  // Need to convert LOCAL x to GLOBAL x in order to correctly calculate
  // PETSC Matrix Index.
  int xoffset = Istart / meshz;
  if( Istart % meshz != 0 )
    throw  BoutException("Petsc index sanity check 3 failed");

  // Calculate the row to be set
  int row_new = x + xshift; // should never be out of range.
  if( !localmesh->firstX() ) row_new += (xoffset - localmesh->xstart);

  // Calculate the column to be set
  int col_new = z + zshift;
  if( col_new < 0 )            col_new += meshz;
  else if( col_new > meshz-1 ) col_new -= meshz;

  // Convert to global indices
  int index = (row_new * meshz) + col_new;

#if CHECK > 2
  if (!finite(ele)) {
    throw BoutException("Non-finite element at x=%d, z=%d, row=%d, col=%d\n",
                        x, z, i, index);
  }
#endif
  
  /* Inserts or adds a block of values into a matrix
   * Input:
   * MatA   - The matrix to set the values in
   * 1      - The number of rows to be set
   * &i     - The global index of the row
   * 1      - The number of columns to be set
   * &index - The global index of the column
   * &ele   - The vlaue to be set
   * INSERT_VALUES replaces existing entries with new values
   */
  MatSetValues(MatA,1,&i,1,&index,&ele,INSERT_VALUES);
}

/*!
 * Set the matrix components of A in Ax=b, solving
 * D*Laplace_perp(x) + (1/C1)Grad_perp(C2)*Grad_perp(x) + Ax = B
 *
 * \note "A" in the equation above is not added here.
 * For calculations of the coefficients, please refer to the user manual.
 *
 * \param[in] x The current x index
 * \param[in] y The current y index
 * \param[in] z The current y index
 * \param[in] coef1  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef2  Convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef3  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef4  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 * \param[in] coef5  Placeholder for convenient variable used to set matrix
 *                   (see manual for details)
 *
 * \param[out] coef1    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef2    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef3    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef4    Convenient variable used to set matrix
 *                      (see manual for details)
 * \param[out] coef5    Convenient variable used to set matrix
 *                      (see manual for details)
 */
void LaplacePetsc::Coeffs( int x, int y, int z, BoutReal &coef1, BoutReal &coef2, BoutReal &coef3, BoutReal &coef4, BoutReal &coef5 ) {

  coef1 = coords->g11(x,y);     // X 2nd derivative coefficient
  coef2 = coords->g33(x,y);     // Z 2nd derivative coefficient
  coef3 = 2.*coords->g13(x,y);  // X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  // If global flag all_terms are set (true by default)
  if (all_terms) {
    coef4 = coords->G1(x,y); // X 1st derivative
    coef5 = coords->G3(x,y); // Z 1st derivative

    ASSERT3(finite(coef4));
    ASSERT3(finite(coef5));
  }

  if(nonuniform) {
    // non-uniform mesh correction
    if((x != 0) && (x != (localmesh->LocalNx-1))) {
      coef4 -= 0.5 * ( ( coords->dx(x+1,y) - coords->dx(x-1,y) ) / SQ(coords->dx(x,y)) ) * coef1; // BOUT-06 term
    }
  }

  if(localmesh->IncIntShear) {
    // d2dz2 term
    coef2 += coords->g11(x,y) * coords->IntShiftTorsion(x,y) * coords->IntShiftTorsion(x,y);
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }

  if (issetD) {
  coef1 *= D(x,y,z);
  coef2 *= D(x,y,z);
  coef3 *= D(x,y,z);
  coef4 *= D(x,y,z);
  coef5 *= D(x,y,z);
  }

  // A second/fourth order derivative term
  if (issetC) {
//   if( (x > 0) && (x < (localmesh->LocalNx-1)) ) //Valid if doing second order derivative, not if fourth: should only be called for xstart<=x<=xend anyway
    if( (x > 1) && (x < (localmesh->LocalNx-2)) ) {
          int zp = z+1;     // z plus 1
          if (zp > meshz-1) zp -= meshz;
          int zm = z-1;     // z minus 1
          if (zm<0) zm += meshz;
          BoutReal ddx_C;
          BoutReal ddz_C;

          if (fourth_order) {
            int zpp = z+2;  // z plus 1 plus 1
            if (zpp > meshz-1) zpp -= meshz;
            int zmm = z-2;  // z minus 1 minus 1
            if (zmm<0) zmm += meshz;
            // Fourth order discretization of C in x
            ddx_C = (-C2(x+2,y,z) + 8.*C2(x+1,y,z) - 8.*C2(x-1,y,z) + C2(x-2,y,z)) / (12.*coords->dx(x,y)*(C1(x,y,z)));
            // Fourth order discretization of C in z
            ddz_C = (-C2(x,y,zpp) + 8.*C2(x,y,zp) - 8.*C2(x,y,zm) + C2(x,y,zmm)) / (12.*coords->dz*(C1(x,y,z)));
          }
          else {
            // Second order discretization of C in x
            ddx_C = (C2(x+1,y,z) - C2(x-1,y,z)) / (2.*coords->dx(x,y)*(C1(x,y,z)));
            // Second order discretization of C in z
            ddz_C = (C2(x,y,zp) - C2(x,y,zm)) / (2.*coords->dz*(C1(x,y,z)));
          }

          coef4 += coords->g11(x,y) * ddx_C + coords->g13(x,y) * ddz_C;
          coef5 += coords->g13(x,y) * ddx_C + coords->g33(x,y) * ddz_C;
        }
    }

  /* Ex and Ez
   * Additional 1st derivative terms to allow for solution field to be
   * components of a vector
   *
   * NB multiply by D or Grad_perp(C)/C as appropriate before passing to
   * setCoefEx()/setCoefEz() because (in principle) both are needed and we
   * don't know how to split them up here
   */
  if (issetE) {
    // These coefficients are 0 by default
    coef4 += Ex(x,y,z);
    coef5 += Ez(x,y,z);
  }
}


void LaplacePetsc::vecToField(Vec xs, FieldPerp &f) {

  ASSERT1(localmesh == f.getMesh());

  f.allocate();
  int i = Istart;
  if(localmesh->firstX())
    {
      for(int x=0; x<localmesh->xstart; x++)
        {
          for(int z=0; z<localmesh->LocalNz; z++)
            {
              PetscScalar val;
              VecGetValues(xs, 1, &i, &val );
              f[x][z] = val;
              i++; // Increment row in Petsc matrix
            }
        }
    }

  for(int x=localmesh->xstart; x <= localmesh->xend; x++)
    {
      for(int z=0; z<localmesh->LocalNz; z++)
        {
          PetscScalar val;
          VecGetValues(xs, 1, &i, &val );
          f[x][z] = val;
          i++; // Increment row in Petsc matrix
        }
    }

  if(localmesh->lastX())
    {
      for(int x=localmesh->xend+1; x<localmesh->LocalNx; x++)
        {
          for(int z=0;z < localmesh->LocalNz; z++)
            {
              PetscScalar val;
              VecGetValues(xs, 1, &i, &val );
              f[x][z] = val;
              i++; // Increment row in Petsc matrix
            }
        }
    }
  ASSERT1(i == Iend);
}

void LaplacePetsc::fieldToVec(const FieldPerp &f, Vec bs) {
  ASSERT1(localmesh == f.getMesh());

  int i = Istart;
  if(localmesh->firstX()) {
    for(int x=0; x<localmesh->xstart; x++) {
      for(int z=0; z<localmesh->LocalNz; z++) {
        PetscScalar val = f[x][z];
        VecSetValues( bs, 1, &i, &val, INSERT_VALUES );
        i++; // Increment row in Petsc matrix
      }
    }
  }

  for(int x=localmesh->xstart; x <= localmesh->xend; x++) {
    for(int z=0; z<localmesh->LocalNz; z++) {
      PetscScalar val = f[x][z];
      VecSetValues( bs, 1, &i, &val, INSERT_VALUES );
      i++; // Increment row in Petsc matrix
    }
  }

  if(localmesh->lastX()) {
    for(int x=localmesh->xend+1; x<localmesh->LocalNx; x++) {
      for(int z=0;z < localmesh->LocalNz; z++) {
        PetscScalar val = f[x][z];
        VecSetValues( bs, 1, &i, &val, INSERT_VALUES );
        i++; // Increment row in Petsc matrix
      }
    }
  }
  ASSERT1(i == Iend);

  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);
}

/// Preconditioner function
int LaplacePetsc::precon(Vec x, Vec y) {
  // Get field to be preconditioned
  FieldPerp xfield;
  vecToField(x, xfield);
  xfield.setIndex(sol.getIndex()); // y index stored in sol variable

  // Call the preconditioner solver
  FieldPerp yfield = pcsolve->solve(xfield);

  // Put result into y
  fieldToVec(yfield, y);
  return 0;
}

#endif // BOUT_HAS_PETSC_3_3
