/*!
 * \file invert_laplace.cxx
 *
 * \brief Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is \f$d*\nabla^2_\perp x + (1./c)\nabla_perp c\cdot\nabla_\perp x + a x = b \f$, where
 * \f$x\f$ and \f$x\f$ are perpendicular (X-Z) or 3D fields,
 * and \f$a\f$ and d are 2D fields. If d is not supplied then it is 1
 *
 * Flags control the boundary conditions (see header file)
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
 */

#include <bout/mesh.hxx>
#include <globals.hxx>
#include <invert_laplace.hxx>
#include <bout_types.hxx>
#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <cmath>
#include <bout/sys/timer.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>
#include <bout/openmpwrap.hxx>

// Implementations:
#include "impls/hypre3d/hypre3d_laplace.hxx"
#include "impls/cyclic/cyclic_laplace.hxx"
#include "impls/iterative_parallel_tri/iterative_parallel_tri.hxx"
#include "impls/multigrid/multigrid_laplace.hxx"
#include "impls/naulin/naulin_laplace.hxx"
#include "impls/pcr/pcr.hxx"
#include "impls/pdd/pdd.hxx"
#include "impls/petsc/petsc_laplace.hxx"
#include "impls/petsc3damg/petsc3damg.hxx"
#include "impls/serial_band/serial_band.hxx"
#include "impls/serial_tri/serial_tri.hxx"
#include "impls/spt/spt.hxx"

/**********************************************************************************
 *                         INITIALISATION AND CREATION
 **********************************************************************************/

/// Laplacian inversion initialisation. Called once at the start to get settings
Laplacian::Laplacian(Options* options, const CELL_LOC loc, Mesh* mesh_in)
    : location(loc), localmesh(mesh_in == nullptr ? bout::globals::mesh : mesh_in) {

  if (options == nullptr) {
    // Use the default options
    options = Options::getRoot()->getSection("laplace");
  }

  output.write("Initialising Laplacian inversion routines\n");

  if (location == CELL_DEFAULT) {
    location = CELL_CENTRE;
  }

  coords = localmesh->getCoordinates(location);

  // Communication option. Controls if asyncronous sends are used
  async_send = (*options)["async"].doc("Use asyncronous MPI send?").withDefault(true);

  const BoutReal filter = (*options)["filter"]
                              .doc("Fraction of Z modes to filter out. Between 0 and 1")
                              .withDefault(0.0);
  const int ncz = localmesh->LocalNz;
  // convert filtering into an integer number of modes
  maxmode = (*options)["maxmode"]
                .doc("The maximum Z mode to solve for")
                .withDefault(ROUND((1.0 - filter) * static_cast<BoutReal>(ncz / 2)));

  // Clamp maxmode between 0 and nz/2
  if (maxmode < 0) {
    maxmode = 0;
  }
  if (maxmode > ncz / 2) {
    maxmode = ncz / 2;
  }

  low_mem = (*options)["low_mem"]
                .doc("If true, reduce the amount of memory used")
                .withDefault(false);

  nonuniform = (*options)["nonuniform"]
                   .doc("Use non-uniform grid corrections? Default is the mesh setting.")
                   .withDefault(coords->non_uniform);

  all_terms =
      (*options)["all_terms"].doc("Include first derivative terms?").withDefault(true);

  global_flags = (*options)["global_flags"].doc("Default flags").withDefault(0);
  inner_boundary_flags = (*options)["inner_boundary_flags"]
                             .doc("Flags to set inner boundary condition")
                             .withDefault(0);
  outer_boundary_flags = (*options)["outer_boundary_flags"]
                             .doc("Flags to set outer boundary condition")
                             .withDefault(0);

  include_yguards = (*options)["include_yguards"]
                        .doc("Solve Laplacian in Y guard cells?")
                        .withDefault(false);

  extra_yguards_lower =
      (*options)["extra_yguards_lower"]
          .doc("Exclude some number of points at the lower boundary, useful for "
               "staggered grids or when boundary conditions make inversion redundant")
          .withDefault(0);
  extra_yguards_upper =
      (*options)["extra_yguards_upper"]
          .doc("Exclude some number of points at the upper boundary, useful for "
               "staggered grids or when boundary conditions make inversion redundant")
          .withDefault(0);
}

std::unique_ptr<Laplacian> Laplacian::instance = nullptr;

Laplacian* Laplacian::defaultInstance() {
  if (instance == nullptr)
    instance = create();
  return instance.get();
}

void Laplacian::cleanup() {
  instance.reset();
}

/**********************************************************************************
 *                                 Solve routines
 **********************************************************************************/

Field3D Laplacian::solve(const Field3D& b) {
  TRACE("Laplacian::solve(Field3D)");

  ASSERT1(b.getLocation() == location);
  ASSERT1(localmesh == b.getMesh());

  Timer timer("invert");
  int ys = localmesh->ystart, ye = localmesh->yend;

  if(localmesh->hasBndryLowerY()) {
    if (include_yguards)
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells

    ys += extra_yguards_lower;
  }
  if(localmesh->hasBndryUpperY()) {
    if (include_yguards)
      ye = localmesh->LocalNy-1; // Contains upper boundary and we are solving in the guard cells

    ye -= extra_yguards_upper;
  }

  Field3D x{emptyFrom(b)};

  int status = 0;
  try {
    for(int jy=ys; jy <= ye; jy++) {
      // 1. Slice b (i.e. take a X-Z plane out of the field)
      // 2. Send it to the solver of the implementation (determined during creation)
      x = solve(sliceXZ(b,jy));
    }
  } catch (const BoutIterationFail&) {
    status = 1;
  }
  BoutParallelThrowRhsFail(status, "Laplacian inversion took too many iterations.");

  return x;
}

// NB: Really inefficient, but functional
Field2D Laplacian::solve(const Field2D& b) {

  ASSERT1(b.getLocation() == location);

  Field3D f = b;
  f = solve(f);
  return DC(f);
}

/*!
 * Performs the laplacian inversion y-slice by y-slice
 *
 * \param[in] b     All the y-slices of b_slice, which is the right hand side
 *                  of the equation A*x_slice = b_slice
 * \param[in] x0    All the y-slices of the variable eventually used to set BC
 *
 * \returns x All the y-slices of x_slice in the equation A*x_slice = b_slice
 */
Field3D Laplacian::solve(const Field3D& b, const Field3D& x0) {
  TRACE("Laplacian::solve(Field3D, Field3D)");

  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());

  Timer timer("invert");

  // Setting the start and end range of the y-slices
  int ys = localmesh->ystart, ye = localmesh->yend;
  if(localmesh->hasBndryLowerY() && include_yguards)
    ys = 0; // Mesh contains a lower boundary
  if(localmesh->hasBndryUpperY() && include_yguards)
    ye = localmesh->LocalNy-1; // Contains upper boundary

  Field3D x{emptyFrom(b)};

  int status = 0;
  try {
    for(int jy=ys; jy <= ye; jy++) {
      // 1. Slice b and x (i.e. take a X-Z plane out of the field)
      // 2. Send them to the solver of the implementation (determined during creation)
      x = solve(sliceXZ(b,jy), sliceXZ(x0,jy));
    }
  } catch (const BoutIterationFail&) {
    status = 1;
  }
  BoutParallelThrowRhsFail(status, "Laplacian inversion took too many iterations.");

  return x; // Return the result of the inversion
}

Field2D Laplacian::solve(const Field2D& b, const Field2D& x0) {
  Field3D f = b, g = x0;
  f = solve(f, g);
  return DC(f);
}

/**********************************************************************************
 *                              MATRIX ELEMENTS
 **********************************************************************************/

void Laplacian::tridagCoefs(int jx, int jy, int jz,
                            dcomplex &a, dcomplex &b, dcomplex &c,
                            const Field2D *ccoef, const Field2D *d,
                            CELL_LOC loc) {

  if (loc == CELL_DEFAULT) loc = location;

  ASSERT1(ccoef == nullptr || ccoef->getLocation() == loc);
  ASSERT1(d == nullptr || d->getLocation() == loc);

  BoutReal kwave=jz*2.0*PI/coords->zlength(); // wave number is 1/[rad]

  tridagCoefs(jx, jy, kwave,
              a, b, c,
              ccoef, d, loc);
}

void Laplacian::tridagCoefs(int jx, int jy, BoutReal kwave,
                            dcomplex &a, dcomplex &b, dcomplex &c,
                            const Field2D *c1coef, const Field2D *c2coef,
                            const Field2D *d, CELL_LOC loc) {
  /* Function: Laplacian::tridagCoef
   * Purpose:  - Set the matrix components of A in Ax=b, solving
   *
   *             D*Laplace_perp(x) + (1/C1)Grad_perp(C2)*Grad_perp(x) + Ax = B
   *
   *             for each fourier component.
   *             NOTE: A in the equation above is not added here.
   *             For calculations of the coefficients, please refer to the user
   *             manual
   *
   * Input:
   * jx        - The current x index
   * jy        - The current y index
   * kwave     - The mode number multiplied with (2*pi)/localmesh->zlength(), where
   *             zlength() is the length of the full z domain (usually 2*pi)
   * a         - Lower diagonal of the tridiagonal matrix. DO NOT CONFUSE WITH A
   * b         - The main diagonal
   * c         - The upper diagonal. DO NOT CONFUSE WITH C (called ccoef here)
   * c1coef    - C1 in the equation above. DO NOT CONFUSE WITH c
   * c2coef    - C2 in the equation above. DO NOT CONFUSE WITH c
   * d         - D in the equation above
   *
   * Output:
   * a         - Lower diagonal of the tridiagonal matrix. DO NOT CONFUSE WITH A
   * b         - The main diagonal
   * c         - The upper diagonal. DO NOT CONFUSE WITH C1, C2 (called c1coef, c2coef
   *             here)
   */

  Coordinates* localcoords;
  if (loc == CELL_DEFAULT) {
    loc = location;
    localcoords = coords;
  } else {
    localcoords = localmesh->getCoordinates(loc);
  }

  BoutReal coef1, coef2, coef3, coef4, coef5;

  coef1=localcoords->g11(jx,jy);     ///< X 2nd derivative coefficient
  coef2=localcoords->g33(jx,jy);     ///< Z 2nd derivative coefficient
  coef3=2.*localcoords->g13(jx,jy);  ///< X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  // If global flag all_terms are set (true by default)
  if(all_terms) {
    coef4 = localcoords->G1(jx,jy); // X 1st derivative
    coef5 = localcoords->G3(jx,jy); // Z 1st derivative
  }

  if (d != nullptr) {
    // Multiply Delp2 component by a factor
    coef1 *= (*d)(jx,jy);
    coef2 *= (*d)(jx,jy);
    coef3 *= (*d)(jx,jy);
    coef4 *= (*d)(jx,jy);
    coef5 *= (*d)(jx,jy);
  }

  if(nonuniform) {
    // non-uniform mesh correction
    if((jx != 0) && (jx != (localmesh->LocalNx-1))) {
      coef4 -= 0.5*((localcoords->dx(jx+1,jy) - localcoords->dx(jx-1,jy))/SQ(localcoords->dx(jx,jy)))*coef1;
    }
  }

  if (c1coef != nullptr) {
    // First derivative terms
    if((jx > 0) && (jx < (localmesh->LocalNx-1))) {
      BoutReal dc2dx_over_c1 = ((*c2coef)(jx+1,jy) - (*c2coef)(jx-1,jy)) / (2.*localcoords->dx(jx,jy)*((*c1coef)(jx,jy)));
      coef4 += localcoords->g11(jx,jy) * dc2dx_over_c1;
      coef5 += localcoords->g13(jx,jy) * dc2dx_over_c1;
    }
  }

  if(localmesh->IncIntShear) {
    // d2dz2 term
    coef2 += localcoords->g11(jx,jy) * localcoords->IntShiftTorsion(jx,jy) * localcoords->IntShiftTorsion(jx,jy);
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }

  coef1 /= SQ(localcoords->dx(jx,jy));
  coef3 /= 2.*localcoords->dx(jx,jy);
  coef4 /= 2.*localcoords->dx(jx,jy);

  a = dcomplex(coef1 - coef4,-kwave*coef3);
  b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,kwave*coef5);
  c = dcomplex(coef1 + coef4,kwave*coef3);
}

/*!
 * Set the matrix components of A in Ax=b
 *
 * This function will
 *      1. Calling tridagCoef, solving
 *
 *         D*Laplace_perp(x) + (1/C1)Grad_perp(C2)*Grad_perp(x) + Ax = B
 *
 *         for each fourier component
 *      2. Set the boundary conditions by setting the first and last rows
 *         properly
 *
 * \param[in] avec      Lower diagonal of the tridiagonal matrix.
 *                      DO NOT CONFUSE WITH "A"
 * \param[in] bvec      The main diagonal
 * \param[in] cvec      The upper diagonal.
 *                      DO NOT CONFUSE WITH "C" (called ccoef here)
 * \param[in] bk        The b in Ax = b
 * \param[in] jy        Index of the current y-slice
 * \param[in] kz        The mode number index
 * \param[in] kwave     The mode number (different from kz only if we are
 *                      taking a part of the z-domain [and not from 0 to 2*pi])
 * \param[in] global_flags          Global flags of the inversion
 * \param[in] inner_boundary_flags  Flags used to set the inner boundary
 * \param[in] outer_boundary_flags  Flags used to set the outer boundary
 * \param[in] a         A in the equation above. DO NOT CONFUSE WITH avec
 * \param[in] c1coef    C1 in the equation above. DO NOT CONFUSE WITH cvec
 * \param[in] c2coef    C2 in the equation above. DO NOT CONFUSE WITH cvec
 * \param[in] d         D in the equation above
 * \param[in] includeguards Whether or not the guard points in x should be used
 *
 * \param[out] avec     Lower diagonal of the tridiagonal matrix.
 *                      DO NOT CONFUSE WITH "A"
 * \param[out] bvec     The main diagonal
 * \param[out] cvec     The upper diagonal.
 *                      DO NOT CONFUSE WITH "C" (called ccoef here)
 */
void Laplacian::tridagMatrix(dcomplex *avec, dcomplex *bvec, dcomplex *cvec,
                             dcomplex *bk, int jy, int kz, BoutReal kwave,
                             int global_flags, int inner_boundary_flags, int outer_boundary_flags,
                             const Field2D *a, const Field2D *c1coef, const Field2D *c2coef,
                             const Field2D *d,
                             bool includeguards) {

  // Better have either both or neither C coefficients
  ASSERT3((c1coef == nullptr and c2coef == nullptr)
          or (c1coef != nullptr and c2coef != nullptr))

  int xs = 0;            // xstart set to the start of x on this processor (including ghost points)
  int xe = localmesh->LocalNx-1;  // xend set to the end of x on this processor (including ghost points)

  // Do not want boundary cells if x is periodic for cyclic solver. Only other solver which
  // works with periodicX is serial_tri, which uses includeguards==true, so the below isn't called.
  if(!includeguards) {
    if(!localmesh->firstX() || localmesh->periodicX)
      xs = localmesh->xstart; // Inner edge is a guard cell
    if(!localmesh->lastX() || localmesh->periodicX)
      xe = localmesh->xend; // Outer edge is a guard cell
  }

  int ncx = xe - xs; // Total number of points in x to be used

  // Setting the width of the boundary.
  // NOTE: The default is a width of (localmesh->xstart) guard cells
  int inbndry = localmesh->xstart, outbndry=localmesh->xstart;

  // If the flags to assign that only one guard cell should be used is set
  if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))  {
    inbndry = outbndry = 1;
  }
  if(inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if(outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  // Loop through our specified x-domain.
  // The boundaries will be set according to the if-statements below.
  for(int ix=0;ix<=ncx;ix++) {
    // Actually set the metric coefficients
    tridagCoefs(xs+ix, jy, kwave, avec[ix], bvec[ix], cvec[ix], c1coef, c2coef, d);
    if (a != nullptr)
      // Add A to bvec (the main diagonal in the matrix)
      bvec[ix] += (*a)(xs+ix,jy);
  }

  // Set the boundary conditions if x is not periodic
  if(!localmesh->periodicX) {
    if(localmesh->firstX()) {
      // INNER BOUNDARY ON THIS PROCESSOR

      // If no user specified value is set on inner boundary, set the first
      // element in b (in the equation AX=b) to 0
      if(!(inner_boundary_flags & (INVERT_RHS | INVERT_SET))) {
        for(int ix=0;ix<inbndry;ix++)
          bk[ix] = 0.;
      }

      // DC i.e. kz = 0 (the offset mode)
      if(kz == 0) {

        if(inner_boundary_flags & INVERT_DC_GRAD && (inner_boundary_flags & INVERT_SET || inner_boundary_flags & INVERT_RHS)) {
          // Zero gradient at inner boundary
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] =  0.;
            bvec[ix] = -1./sqrt(coords->g_11(ix,jy))/coords->dx(ix,jy);
            cvec[ix] =  1./sqrt(coords->g_11(ix,jy))/coords->dx(ix,jy);
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRAD) {
          // Zero gradient at inner boundary
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] =  0.;
            bvec[ix] = -1.;
            cvec[ix] =  1.;
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRADPAR) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  1.0/sqrt(coords->g_22(ix,jy));
            cvec[ix] = -1.0/sqrt(coords->g_22(ix+1,jy));
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRADPARINV) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  sqrt(coords->g_22(ix,jy));
            cvec[ix] = -sqrt(coords->g_22(ix+1,jy));
          }
        }
        else if (inner_boundary_flags & INVERT_DC_LAP) {
          // Decaying boundary conditions
          BoutReal k = 0.0;
          if (a != nullptr) {
            BoutReal ksq = -((*a)(inbndry, jy));
            if(ksq < 0.0)
              throw BoutException("ksq must be positive");
            k = sqrt(ksq);
          }
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] =  0.;
            bvec[ix] =  1.;
            cvec[ix] = -exp(-k*coords->dx(ix,jy)/sqrt(coords->g11(ix,jy)));
          }
        }
        else if (inner_boundary_flags & INVERT_IN_CYLINDER){
          // Condition for inner radial boundary for cylindrical coordinates
          /* Explanation:
           * The discrete fourier transform is defined as
           *
           * F(x,y)_k = (1/N) sum_{Z=0}^{N-1} f(x,y)_Z exp(-2pi i k Z/N)
           *
           * where
           *
           * k = mode number
           * N = total number of points in the periodic direction
           * Z = index number for point in Z
           * i = imaginary number
           *
           * If the boundary is located at rho=0, and if the ghost point is
           * chosen half between grid cells, then (as there is no real boundary
           * in the center of a cylinder), the ghost point of a innermost inner
           * point in rho will have its ghost point as the innermost inner
           * point diamatrically opposite (i.e. the ghost point is the current
           * point with a rotation of pi [rotation of N/2 in Z if we use index
           * values]). The fourier components at the ghost point will thus be
           * phase shifted with respect to the innermost point under
           * consideration, and for the ghost point we can write
           *
           * F(x,y)_k_ghost
           * = (1/N) sum_{Z=0}^{N-1} f(x,y)_(Z+(N/2)) exp(-2pi i k (Z + N/2)/N)
           * = (1/N) sum_{Z=0}^{N-1} f(x,y)_(Z+(N/2)) exp(-2pi i k Z/N)
           *                                          exp(-pi i k)
           * = F(x,y)_k exp(-pi i k)
           *
           * so
           *
           * F(x,y)_k_ghost = F(x,y)_k exp(-pi i k)
           *                =   F(x,y)_k     if k is even
           *                = - F(x,y)_k     if k is odd
           */
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] = 0.;
            bvec[ix] = 1. ;
            cvec[ix] = -1. ;
          }
          // Zero value at inner boundary
        }
        else {
          // Order 2 dirichlet BC (boundary half between points)
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] = 0.;
            bvec[ix] = 0.5;
            cvec[ix] = 0.5;
          }
        }
      }
      // AC i.e. kz =/= 0 (all other modes than the offset mode)
      else {

        if(inner_boundary_flags & INVERT_AC_GRAD && (inner_boundary_flags & INVERT_SET || inner_boundary_flags & INVERT_RHS)) {
          // Zero gradient at inner boundary
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] = dcomplex(0.,0.);
            bvec[ix] = dcomplex(-1.,0.)/sqrt(coords->g_11(ix,jy))/coords->dx(ix,jy);
            cvec[ix] = dcomplex(1.,0.)/sqrt(coords->g_11(ix,jy))/coords->dx(ix,jy);
          }
        }
        else if(inner_boundary_flags & INVERT_AC_GRAD) {
          // Zero gradient at inner boundary
          for (int ix=0;ix<inbndry;ix++){
            avec[ix]=dcomplex(0.,0.);
            bvec[ix]=dcomplex(-1.,0.);
            cvec[ix]=dcomplex(1.,0.);
          }
        }
        else if(inner_boundary_flags & INVERT_AC_LAP) {
          // Use decaying zero-Laplacian solution in the boundary
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] = 0.0;
            bvec[ix] = 1.0;
            cvec[ix] = -exp(-1.0*sqrt(coords->g33(ix,jy)/coords->g11(ix,jy))*kwave*coords->dx(ix,jy));
          }
        }
        else if (inner_boundary_flags & INVERT_IN_CYLINDER) {
          // Condition for inner radial boundary for cylindrical coordinates
          // Explanation under "if (inner_boundary_flags & INVERT_IN_CYLINDER)"
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] = 0.;
            bvec[ix] = 1.;

            if ((kz % 2) == 0){
              cvec[ix] = -1.;
            }
            else {
              cvec[ix] = 1.;
            }
          }
        }
        else {
          // Order 2 dirichlet BC (boundary half between points)
          // Zero value at inner boundary or INVERT_IN_SET
          for (int ix=0;ix<inbndry;ix++){
            avec[ix]=dcomplex(0.,0.);
            bvec[ix]=dcomplex(0.5,0.);
            cvec[ix]=dcomplex(0.5,0.);
          }
        }
      }
    }
    if(localmesh->lastX()) {
      // OUTER BOUNDARY ON THIS PROCESSOR

      // If no user specified value is set on outer boundary, set the last
      // element in b (in the equation AX=b) to 0
      if(!(outer_boundary_flags & (INVERT_RHS | INVERT_SET))) {
        for (int ix=0;ix<outbndry;ix++) {
          bk[ncx-ix] = 0.;
        }
      }

      // DC i.e. kz = 0 (the offset mode)
      if(kz==0) {

        if(outer_boundary_flags & INVERT_DC_GRAD && ( outer_boundary_flags & INVERT_SET || outer_boundary_flags & INVERT_RHS)) {
          // Zero gradient at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            avec[ncx-ix]=dcomplex(-1.,0.)/sqrt(coords->g_11(ncx-ix,jy))/coords->dx(ncx-ix,jy);
            bvec[ncx-ix]=dcomplex(1.,0.)/sqrt(coords->g_11(ncx-ix,jy))/coords->dx(ncx-ix,jy);
            cvec[ncx-ix]=dcomplex(0.,0.);
          }
        }
        else if(outer_boundary_flags & INVERT_DC_GRAD) {
          // Zero gradient at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            avec[ncx-ix]=dcomplex(1.,0.);
            bvec[ncx-ix]=dcomplex(-1.,0.);
            cvec[ncx-ix]=dcomplex(0.,0.);
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRADPAR) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ncx-ix] =  1.0/sqrt(coords->g_22(ncx-ix+1,jy));
            bvec[ncx-ix] = -1.0/sqrt(coords->g_22(ncx-ix,jy));
            cvec[ncx-ix] =  0.0;
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRADPARINV) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ncx-ix] =  sqrt(coords->g_22(ncx-ix-1,jy));
            bvec[ncx-ix] = -sqrt(coords->g_22(ncx-ix,jy));
            cvec[ncx-ix] =  0.0;
          }
        }
        else if (inner_boundary_flags & INVERT_DC_LAP) {
          // Decaying boundary conditions
          BoutReal k = 0.0;
          if (a != nullptr) {
            BoutReal ksq = -((*a)(inbndry, jy));
            if(ksq < 0.0)
              throw BoutException("ksq must be positive");
            k = sqrt(ksq);
          }
          for (int ix=0;ix<inbndry;ix++){
            cvec[ncx-ix] =  0.;
            bvec[ncx-ix] =  1.;
            avec[ncx-ix] = -exp(-k*coords->dx(ncx-ix,jy)/sqrt(coords->g11(ncx-ix,jy)));
          }
        }
        else {
          // Order 2 dirichlet BC (boundary half between points)
          // Zero value at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            cvec[ncx-ix]=dcomplex(0.,0.);
            bvec[ncx-ix]=dcomplex(0.5,0.);
            avec[ncx-ix]=dcomplex(0.5,0.);
          }
        }
      }
      // AC i.e. kz =/= 0 (all other modes than the offset mode)
      else {

        if(outer_boundary_flags & INVERT_AC_GRAD && ( outer_boundary_flags & INVERT_SET || outer_boundary_flags & INVERT_RHS)) {
          // Zero gradient at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            avec[ncx-ix]=dcomplex(-1.,0.)/sqrt(coords->g_11(ncx-ix,jy))/coords->dx(ncx-ix,jy);
            bvec[ncx-ix]=dcomplex(1.,0.)/sqrt(coords->g_11(ncx-ix,jy))/coords->dx(ncx-ix,jy);
            cvec[ncx-ix]=dcomplex(0.,0.);
          }
        }
        else if(outer_boundary_flags & INVERT_AC_GRAD) {
          // Zero gradient at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            avec[ncx-ix]=dcomplex(1.,0.);
            bvec[ncx-ix]=dcomplex(-1.,0.);
            cvec[ncx-ix]=dcomplex(0.,0.);
          }
        }
        else if(outer_boundary_flags & INVERT_AC_LAP) {
          // Use decaying zero-Laplacian solution in the boundary
          for (int ix=0;ix<outbndry;ix++) {
            avec[ncx-ix] = -exp(-1.0*sqrt(coords->g33(xe-ix,jy)/coords->g11(xe-ix,jy))*kwave*coords->dx(xe-ix,jy));
            bvec[ncx-ix] = 1.0;
            cvec[ncx-ix] = 0.0;
          }
        }
        else {
          // Order 2 dirichlet BC (boundary half between points)
          // Zero value at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            cvec[ncx-ix]=dcomplex(0.,0.);
            bvec[ncx-ix]=dcomplex(0.5,0.);
            avec[ncx-ix]=dcomplex(0.5,0.);
          }
        }
      }
    }
  }
}

/**********************************************************************************
 *                              LEGACY INTERFACE
 *
 * These functions are depreciated, and will be removed in future
 **********************************************************************************/

/// Returns the coefficients for a tridiagonal matrix for laplace. Used by Delp2 too
void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c,
                          const Field2D *ccoef, const Field2D *d, CELL_LOC loc) {
  Laplacian::defaultInstance()->tridagCoefs(jx,jy, jz, a, b, c, ccoef, d, loc);
}

constexpr decltype(LaplaceFactory::type_name) LaplaceFactory::type_name;
constexpr decltype(LaplaceFactory::section_name) LaplaceFactory::section_name;
constexpr decltype(LaplaceFactory::option_name) LaplaceFactory::option_name;
constexpr decltype(LaplaceFactory::default_type) LaplaceFactory::default_type;
