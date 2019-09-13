/**************************************************************************
 * 3D Laplacian Solver
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

#include "petsc3damg.hxx"

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

LaplacePetsc3dAmg::LaplacePetsc3dAmg(Options *opt, const CELL_LOC loc, Mesh *mesh_in) :
  Laplacian(opt, loc, mesh_in),
  A(0.0), C1(1.0), C2(1.0), D(1.0), Ex(0.0), Ez(0.0),
  issetD(false), issetC(false), issetE(false)
{

  // Provide basic initialisation of field coefficients, etc.
  // Get relevent options from user input
  // Initialise PETSc objects
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
      throw BoutException("LaplacePetsc3dAmg does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
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
    throw BoutException("Error in MPI_Allreduce during LaplacePetsc3dAmg initialisation");

  // Calculate total (physical) grid dimensions
  meshz = localmesh->LocalNz;
  meshx = size / meshz;

  // Get 4th order solver switch
  opts->get("fourth_order", fourth_order, false);

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
}

const Field3D LaplacePetsc3dAmg::solve(const Field3D &b_in, const Field3D &x0) {
  // Copy initial guess and RHS into PETSc vectors
  // Assemble matrix operator
    // Handle internal values (easy)
    // Handle boundary conditions (which ones?)
  // Configure preconditioner (or do this in constructor?)
}

const Field2D LaplacePetsc3dAmg::solve(const Field2D &b) {
  return Laplacian::solve(b);
}

PetscMatrix<Field3D>& LaplacePetsc3dAmg::getMatrix3D() {
}

PetscMatrix<Field2D>& LaplacePetsc3dAmg::getMatrix2D() {
}


#endif // BOUT_HAS_PETSC_3_3
