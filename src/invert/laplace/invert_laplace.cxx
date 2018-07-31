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

#include "laplacefactory.hxx"

/**********************************************************************************
 *                         INITIALISATION AND CREATION
 **********************************************************************************/

/// Laplacian inversion initialisation. Called once at the start to get settings
Laplacian::Laplacian(Options *options) {

  if (options == nullptr) {
    // Use the default options
    options = Options::getRoot()->getSection("laplace");
  }

  output.write("Initialising Laplacian inversion routines\n");

  // Communication option. Controls if asyncronous sends are used
  options->get("async", async_send, true);

  BoutReal filter; ///< Fraction of Z modes to filter out. Between 0 and 1
  OPTION(options, filter, 0.0);
  int ncz = mesh->LocalNz;
  // convert filtering into an integer number of modes
  maxmode = ROUND((1.0 - filter) * static_cast<BoutReal>(ncz / 2));
  // Can be overriden by max_mode option
  OPTION(options, maxmode, maxmode);
  if(maxmode < 0) maxmode = 0;
  if(maxmode > ncz/2) maxmode = ncz/2;

  OPTION(options, low_mem, false);

  OPTION(options, nonuniform,
         mesh->coordinates()->non_uniform); // Default is the mesh setting

  OPTION(options, all_terms, true); // Include first derivative terms

  if (options->isSet("flags")) {
    if ( options->isSet("global_flags") || options->isSet("inner_boundary_flags") || options->isSet("outer_boundary_flags") ) {
      throw BoutException("Should not use old flags as well as new global_flags/inner_boundary_flags/outer_boundary_flags");
    }
    int flags;
    OPTION(options, flags, 0);
    setFlags(flags);
  }
  else {
    OPTION(options, global_flags, 0);
    OPTION(options, inner_boundary_flags, 0);
    OPTION(options, outer_boundary_flags, 0);
  }

  OPTION(options, include_yguards, false);

  OPTION2(options, extra_yguards_lower, extra_yguards_upper, 0);
}

Laplacian* Laplacian::create(Options *opts) {
  // Factory pattern:
  // 1. getInstance() is making an instance of LaplacianFactory
  // 2. createLaplacian() is accessing this instance and returning a Laplacian
  //    form one of the child classes of the Laplacian (the laplace solver
  //    implementations)
  return LaplaceFactory::getInstance()->createLaplacian(opts);
}

Laplacian *Laplacian::instance = nullptr;

Laplacian* Laplacian::defaultInstance() {
  if (instance == nullptr)
    instance = create();
  return instance;
}

void Laplacian::cleanup() {
  if (instance == nullptr)
    return;
  delete instance;
  instance = nullptr;
}

/**********************************************************************************
 *                                 Solve routines
 **********************************************************************************/

const Field3D Laplacian::solve(const Field3D &b) {
  TRACE("Laplacian::solve(Field3D)");
  Mesh *mesh = b.getMesh();

  Timer timer("invert");
  int ys = mesh->ystart, ye = mesh->yend;

  if(mesh->hasBndryLowerY()) {
    if (include_yguards)
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells

    ys += extra_yguards_lower;
  }
  if(mesh->hasBndryUpperY()) {
    if (include_yguards)
      ye = mesh->LocalNy-1; // Contains upper boundary and we are solving in the guard cells

    ye -= extra_yguards_upper;
  }

  Field3D x(mesh);
  x.allocate();

  int status = 0;
  try {
    for(int jy=ys; jy <= ye; jy++) {
      // 1. Slice b (i.e. take a X-Z plane out of the field)
      // 2. Send it to the solver of the implementation (determined during creation)
      x = solve(sliceXZ(b,jy));
    }
  } catch (BoutIterationFail &itfail) {
    status = 1;
  }
  BoutParallelThrowRhsFail(status, "Laplacian inversion took too many iterations.");

  x.setLocation(b.getLocation());

  return x;
}

// NB: Really inefficient, but functional
const Field2D Laplacian::solve(const Field2D &b) {
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
const Field3D Laplacian::solve(const Field3D &b, const Field3D &x0) {
  TRACE("Laplacian::solve(Field3D, Field3D)");

  Timer timer("invert");

  Mesh *mesh = b.getMesh();
  // Setting the start and end range of the y-slices
  int ys = mesh->ystart, ye = mesh->yend;
  if(mesh->hasBndryLowerY() && include_yguards)
    ys = 0; // Mesh contains a lower boundary
  if(mesh->hasBndryUpperY() && include_yguards)
    ye = mesh->LocalNy-1; // Contains upper boundary

  Field3D x(mesh);
  x.allocate();

  int status = 0;
  try {
    for(int jy=ys; jy <= ye; jy++) {
      // 1. Slice b and x (i.e. take a X-Z plane out of the field)
      // 2. Send them to the solver of the implementation (determined during creation)
      x = solve(sliceXZ(b,jy), sliceXZ(x0,jy));
    }
  } catch (BoutIterationFail &itfail) {
    status = 1;
  }
  BoutParallelThrowRhsFail(status, "Laplacian inversion took too many iterations.");

  x.setLocation(b.getLocation());

  return x; // Return the result of the inversion
}

const Field2D Laplacian::solve(const Field2D &b, const Field2D &x0) {
  Field3D f = b, g = x0;
  f = solve(f, g);
  return DC(f);
}

/**********************************************************************************
 *                              MATRIX ELEMENTS
 **********************************************************************************/

void Laplacian::tridagCoefs(int jx, int jy, int jz,
                            dcomplex &a, dcomplex &b, dcomplex &c,
                            const Field2D *ccoef, const Field2D *d, CELL_LOC outloc) {

  Coordinates *coord = mesh->coordinates();

  BoutReal kwave=jz*2.0*PI/coord->zlength(); // wave number is 1/[rad]

  tridagCoefs(jx, jy, kwave,
              a, b, c,
              ccoef, d, outloc);
}

void Laplacian::tridagCoefs(int jx, int jy, BoutReal kwave,
                            dcomplex &a, dcomplex &b, dcomplex &c,
                            const Field2D *ccoef, const Field2D *d, CELL_LOC outloc) {
  /* Function: Laplacian::tridagCoef
   * Purpose:  - Set the matrix components of A in Ax=b, solving
   *
   *             D*Laplace_perp(x) + (1/C)Grad_perp(C)*Grad_perp(x) + Ax = B
   *
   *             for each fourier component.
   *             NOTE: A in the equation above is not added here.
   *             For calculations of the coefficients, please refer to the user
   *             manual
   *
   * Input:
   * jx        - The current x index
   * jy        - The current y index
   * kwave     - The mode number multiplied with (2*pi)/mesh->zlength(), where
   *             zlength() is the length of the full z domain (usually 2*pi)
   * a         - Lower diagonal of the tridiagonal matrix. DO NOT CONFUSE WITH A
   * b         - The main diagonal
   * c         - The upper diagonal. DO NOT CONFUSE WITH C (called ccoef here)
   * ccoef     - C in the equation above. DO NOT CONFUSE WITH c
   * d         - D in the equation above
   *
   * Output:
   * a         - Lower diagonal of the tridiagonal matrix. DO NOT CONFUSE WITH A
   * b         - The main diagonal
   * c         - The upper diagonal. DO NOT CONFUSE WITH C (called ccoef here)
   */
  BoutReal coef1, coef2, coef3, coef4, coef5;

  Coordinates *coord = mesh->coordinates();

  if (outloc == CELL_DEFAULT) {
    outloc = CELL_CENTRE;
  }
  const Field2D & g11 = coord->g11.get(outloc);
  const Field2D & g33 = coord->g33.get(outloc);
  const Field2D & g13 = coord->g13.get(outloc);
  const Field2D & G1 = coord->G1.get(outloc);
  const Field2D & G3 = coord->G3.get(outloc);
  const Field2D & dx = coord->dx.get(outloc);

  coef1=g11(jx,jy);     ///< X 2nd derivative coefficient
  coef2=g33(jx,jy);     ///< Z 2nd derivative coefficient
  coef3=2.*g13(jx,jy);  ///< X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  // If global flag all_terms are set (true by default)
  if(all_terms) {
    coef4 = G1(jx,jy); // X 1st derivative
    coef5 = G3(jx,jy); // Z 1st derivative
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
    if((jx != 0) && (jx != (mesh->LocalNx-1))) {
      coef4 -= 0.5*((dx(jx+1,jy) - dx(jx-1,jy))/SQ(dx(jx,jy)))*coef1;
    }
  }

  if (ccoef != nullptr) {
    // A first order derivative term
    if((jx > 0) && (jx < (mesh->LocalNx-1)))
      coef4 += g11(jx,jy) * ((*ccoef)(jx+1,jy) - (*ccoef)(jx-1,jy)) / (2.*dx(jx,jy)*((*ccoef)(jx,jy)));
  }

  if(mesh->IncIntShear) {
    Field2D IntShiftTorsion = coord->IntShiftTorsion.get(outloc);
    // d2dz2 term
    coef2 += g11(jx,jy) * IntShiftTorsion(jx,jy) * IntShiftTorsion(jx,jy);
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }

  coef1 /= SQ(dx(jx,jy));
  coef3 /= 2.*dx(jx,jy);
  coef4 /= 2.*dx(jx,jy);

  a = dcomplex(coef1 - coef4,-kwave*coef3);
  b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,kwave*coef5);
  c = dcomplex(coef1 + coef4,kwave*coef3);
}

/// Sets the coefficients for parallel tridiagonal matrix inversion
/*!
 * Uses the laplace_tridag_coefs routine above to fill a matrix [kz][ix] of coefficients
 */
void Laplacian::tridagMatrix(dcomplex **avec, dcomplex **bvec, dcomplex **cvec,
                             dcomplex **bk, int jy, int global_flags, int inner_boundary_flags, int outer_boundary_flags,
                             const Field2D *a, const Field2D *ccoef,
                             const Field2D *d) {

  Coordinates *coord = mesh->coordinates();

  BOUT_OMP(parallel for)
  for(int kz = 0; kz <= maxmode; kz++) {
    BoutReal kwave=kz*2.0*PI/coord->zlength(); // wave number is 1/[rad]

    tridagMatrix(avec[kz], bvec[kz], cvec[kz],
                 bk[kz],
                 jy,
                 kz, kwave,
                 global_flags, inner_boundary_flags, outer_boundary_flags,
                 a, ccoef, d);
  }
}

/*!
 * Set the matrix components of A in Ax=b
 *
 * This function will
 *      1. Calling tridagCoef, solving
 *
 *         D*Laplace_perp(x) + (1/C)Grad_perp(C)*Grad_perp(x) + Ax = B
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
 * \param[in] ccoef     C in the equation above. DO NOT CONFUSE WITH cvec
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
                             const Field2D *a, const Field2D *ccoef,
                             const Field2D *d,
                             bool includeguards) {
  int xs = 0;            // xstart set to the start of x on this processor (including ghost points)
  int xe = mesh->LocalNx-1;  // xend set to the end of x on this processor (including ghost points)

  Coordinates *coord = mesh->coordinates();

  // Do not want boundary cells if x is periodic for cyclic solver. Only other solver which
  // works with periodicX is serial_tri, which uses includeguards==true, so the below isn't called.
  if(!includeguards) {
    if(!mesh->firstX() || mesh->periodicX)
      xs = mesh->xstart; // Inner edge is a guard cell
    if(!mesh->lastX() || mesh->periodicX)
      xe = mesh->xend; // Outer edge is a guard cell
  }

  int ncx = xe - xs; // Total number of points in x to be used

  // Setting the width of the boundary.
  // NOTE: The default is a width of (mesh->xstart) guard cells
  int inbndry = mesh->xstart, outbndry=mesh->xstart;

  // If the flags to assign that only one guard cell should be used is set
  if((global_flags & INVERT_BOTH_BNDRY_ONE) || (mesh->xstart < 2))  {
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
    tridagCoefs(xs+ix, jy, kwave, avec[ix], bvec[ix], cvec[ix], ccoef, d);
    if (a != nullptr)
      // Add A to bvec (the main diagonal in the matrix)
      bvec[ix] += (*a)(xs+ix,jy);
  }

  // Set the boundary conditions if x is not periodic
  if(!mesh->periodicX) {
    if(mesh->firstX()) {
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
            bvec[ix] = -1./sqrt(coord->g_11(ix,jy))/coord->dx(ix,jy);
            cvec[ix] =  1./sqrt(coord->g_11(ix,jy))/coord->dx(ix,jy);
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
            bvec[ix] =  1.0/sqrt(coord->g_22(ix,jy));
            cvec[ix] = -1.0/sqrt(coord->g_22(ix+1,jy));
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRADPARINV) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  sqrt(coord->g_22(ix,jy));
            cvec[ix] = -sqrt(coord->g_22(ix+1,jy));
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
            cvec[ix] = -exp(-k*coord->dx(ix,jy)/sqrt(coord->g11(ix,jy)));
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
            bvec[ix] = dcomplex(-1.,0.)/sqrt(coord->g_11(ix,jy))/coord->dx(ix,jy);
            cvec[ix] = dcomplex(1.,0.)/sqrt(coord->g_11(ix,jy))/coord->dx(ix,jy);
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
            cvec[ix] = -exp(-1.0*sqrt(coord->g33(ix,jy)/coord->g11(ix,jy))*kwave*coord->dx(ix,jy));
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
    if(mesh->lastX()) {
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
            avec[ncx-ix]=dcomplex(-1.,0.)/sqrt(coord->g_11(ncx-ix,jy))/coord->dx(ncx-ix,jy);
            bvec[ncx-ix]=dcomplex(1.,0.)/sqrt(coord->g_11(ncx-ix,jy))/coord->dx(ncx-ix,jy);
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
            avec[ncx-ix] =  1.0/sqrt(coord->g_22(ncx-ix+1,jy));
            bvec[ncx-ix] = -1.0/sqrt(coord->g_22(ncx-ix,jy));
            cvec[ncx-ix] =  0.0;
          }
        }
        else if(inner_boundary_flags & INVERT_DC_GRADPARINV) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ncx-ix] =  sqrt(coord->g_22(ncx-ix-1,jy));
            bvec[ncx-ix] = -sqrt(coord->g_22(ncx-ix,jy));
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
            avec[ncx-ix] = -exp(-k*coord->dx(ncx-ix,jy)/sqrt(coord->g11(ncx-ix,jy)));
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
            avec[ncx-ix]=dcomplex(-1.,0.)/sqrt(coord->g_11(ncx-ix,jy))/coord->dx(ncx-ix,jy);
            bvec[ncx-ix]=dcomplex(1.,0.)/sqrt(coord->g_11(ncx-ix,jy))/coord->dx(ncx-ix,jy);
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
            avec[ncx-ix] = -exp(-1.0*sqrt(coord->g33(xe-ix,jy)/coord->g11(xe-ix,jy))*kwave*coord->dx(xe-ix,jy));
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
                          const Field2D *ccoef, const Field2D *d, CELL_LOC outloc) {
  Laplacian::defaultInstance()->tridagCoefs(jx,jy, jz, a, b, c, ccoef, d, outloc);
}

int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {

  Laplacian *lap = Laplacian::defaultInstance();

  if (a != nullptr) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);

  if (c != nullptr) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);

  if (d != nullptr) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);

  lap->setFlags(flags);

  x = lap->solve(b);

  x.setLocation(b.getLocation());

  return 0;
}

int invert_laplace(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {

  Timer timer("invert"); ///< Start timer

  Laplacian *lap = Laplacian::defaultInstance();

  if (a != nullptr) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);

  if (c != nullptr) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);

  if (d != nullptr) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);

  lap->setFlags(flags);

  x.allocate(); // Make sure x is allocated

  x = lap->solve(b, x);

  x.setLocation(b.getLocation());

  return 0;
}
const Field3D invert_laplace(const Field3D &b, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {

  Timer timer("invert"); ///< Start timer

  Laplacian *lap = Laplacian::defaultInstance();

  if (a != nullptr) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);

  if (c != nullptr) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);

  if (d != nullptr) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);

  lap->setFlags(flags);

  Field3D x = lap->solve(b);

  x.setLocation(b.getLocation());

  return x;
}

// setFlags routine for backwards compatibility with old monolithic flags
void Laplacian::setFlags(int flags) {
  global_flags = 0;
  inner_boundary_flags = 0;
  outer_boundary_flags = 0;
  if (flags & 1)
    inner_boundary_flags += INVERT_DC_GRAD;
  if (flags & 2)
    inner_boundary_flags += INVERT_AC_GRAD;
  if (flags & 4)
    outer_boundary_flags += INVERT_DC_GRAD;
  if (flags & 8)
    outer_boundary_flags += INVERT_AC_GRAD;
  if (flags & 16)
    global_flags += INVERT_ZERO_DC;
  if (flags & 32)
    global_flags += INVERT_START_NEW;
  if (flags & 64)
    global_flags += INVERT_BOTH_BNDRY_ONE; // Sets the width of the boundary to 1
  if (flags & 128)
    global_flags += INVERT_4TH_ORDER; // Use band solver for 4th order in x
  if (flags & 256)
    inner_boundary_flags += INVERT_AC_LAP;
  if (flags & 512)
    outer_boundary_flags += INVERT_AC_LAP;
  if (flags & 1024)
    inner_boundary_flags += INVERT_SYM; // Use symmetry to enforce either zero-value or zero-gradient
  if (flags & 2048)
    outer_boundary_flags += INVERT_SYM; // Use symmetry to enforce either zero-value or zero-gradient
  if (flags & 4096)
    inner_boundary_flags += INVERT_SET; // Set inner boundary
  if (flags & 8192)
    outer_boundary_flags += INVERT_SET; // Set outer boundary
  if (flags & 16384)
    inner_boundary_flags += INVERT_RHS; // Use input value in RHS at inner boundary
  if (flags & 32768)
    outer_boundary_flags += INVERT_RHS; // Use input value in RHS at outer boundary
  if (flags & 65536)
    global_flags += INVERT_KX_ZERO; // Zero the kx=0, n = 0 component
  if (flags & 131072)
    inner_boundary_flags += INVERT_DC_LAP;
  if (flags & 262144)
    inner_boundary_flags += INVERT_BNDRY_ONE;
  if (flags & 524288)
    outer_boundary_flags += INVERT_BNDRY_ONE;
  if (flags & 1048576)
    inner_boundary_flags += INVERT_DC_GRADPAR;
  if (flags & 2097152)
    inner_boundary_flags += INVERT_DC_GRADPARINV;
  if (flags & 4194304)
    inner_boundary_flags += INVERT_IN_CYLINDER;
}
