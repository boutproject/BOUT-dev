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

#include "laplacefactory.hxx"

/**********************************************************************************
 *                         INITIALISATION AND CREATION
 **********************************************************************************/

/// Laplacian inversion initialisation. Called once at the start to get settings
Laplacian::Laplacian(Options *options) {

  if(options == NULL) {
    // Use the default options
    options = Options::getRoot()->getSection("laplace");
  }
  
  output.write("Initialising Laplacian inversion routines\n");
  
  // Communication option. Controls if asyncronous sends are used
  options->get("async", async_send, true);
  
  BoutReal filter; ///< Fraction of Z modes to filter out. Between 0 and 1
  OPTION(options, filter, 0.2);
  int ncz = mesh->ngz-1;
  // convert filtering into an integer number of modes
  maxmode = ROUND((1.0 - filter) * ((double) (ncz / 2)));
  // Can be overriden by max_mode option
  OPTION(options, maxmode, maxmode);
  if(maxmode < 0) maxmode = 0;
  if(maxmode > ncz/2) maxmode = ncz/2;
  
  OPTION3(options, low_mem, all_terms, nonuniform, false);
  
  OPTION(options, flags, 0);
  
  OPTION(options, include_yguards, true);
  
  OPTION2(options, extra_yguards_lower, extra_yguards_upper, 0);
}

Laplacian* Laplacian::create(Options *opts) {
  return LaplaceFactory::getInstance()->createLaplacian(opts);
}

Laplacian* Laplacian::instance = NULL;

Laplacian* Laplacian::defaultInstance() {
  if(instance == NULL)
    instance = create();
  return instance;
}

/**********************************************************************************
 *                                 Solve routines
 **********************************************************************************/

const Field3D Laplacian::solve(const Field3D &b) {
  Timer timer("invert");
#ifdef CHECK
  msg_stack.push("Laplacian::solve(Field3D)");
#endif
  int ys = mesh->ystart, ye = mesh->yend;

  if(mesh->hasBndryLowerY()) {
    if (include_yguards)
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells
    
    ys += extra_yguards_lower;
  }
  if(mesh->hasBndryUpperY()) {
    if (include_yguards)
      ye = mesh->ngy-1; // Contains upper boundary and we are solving in the guard cells
      
    ye -= extra_yguards_upper;
  }

  Field3D x;
  x.allocate();
  
  int status = 0;
  try {
    for(int jy=ys; jy <= ye; jy++) {
      x = solve(b.slice(jy));
    }
  }
  catch (BoutIterationFail itfail) {
    status = 1;
  }
  BoutParallelThrowRhsFail(status, "Laplacian inversion took too many iterations.");
  
#ifdef CHECK
  msg_stack.pop();
#endif

  x.setLocation(b.getLocation());

  return x;
}

// NB: Really inefficient, but functional
const Field2D Laplacian::solve(const Field2D &b) {
  Field3D f = b;
  f = solve(f);
  return f.DC();
}

const Field3D Laplacian::solve(const Field3D &b, const Field3D &x0) {
  Timer timer("invert");
#ifdef CHECK
  msg_stack.push("Laplacian::solve(Field3D, Field3D)");
#endif

  int ys = mesh->ystart, ye = mesh->yend;
  if(mesh->hasBndryLowerY() && include_yguards)
    ys = 0; // Mesh contains a lower boundary
  if(mesh->hasBndryUpperY() && include_yguards)
    ye = mesh->ngy-1; // Contains upper boundary
  
  Field3D x;
  x.allocate();
  
  int status = 0;
  try {
    for(int jy=ys; jy <= ye; jy++) {
      x = solve(b.slice(jy), x0.slice(jy));
    }
  }
  catch (BoutIterationFail itfail) {
    status = 1;
  }
  BoutParallelThrowRhsFail(status, "Laplacian inversion took too many iterations.");
  
#ifdef CHECK
  msg_stack.pop();
#endif

  x.setLocation(b.getLocation());

  return x;
}

const Field2D Laplacian::solve(const Field2D &b, const Field2D &x0) {
  Field3D f = b, g = x0;
  f = solve(f, g);
  return f.DC();
}

/**********************************************************************************
 *                              MATRIX ELEMENTS
 **********************************************************************************/

void Laplacian::tridagCoefs(int jx, int jy, int jz, 
                            dcomplex &a, dcomplex &b, dcomplex &c, 
                            const Field2D *ccoef, const Field2D *d) {
  
  BoutReal kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
  
  tridagCoefs(jx, jy, kwave, 
              a, b, c, 
              ccoef, d);
}

void Laplacian::tridagCoefs(int jx, int jy, BoutReal kwave, 
                            dcomplex &a, dcomplex &b, dcomplex &c, 
                            const Field2D *ccoef, const Field2D *d) {
  
  BoutReal coef1, coef2, coef3, coef4, coef5;
  
  coef1=mesh->g11[jx][jy];     ///< X 2nd derivative coefficient
  coef2=mesh->g33[jx][jy];     ///< Z 2nd derivative coefficient
  coef3=2.*mesh->g13[jx][jy];  ///< X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  if(all_terms) {
    coef4 = mesh->G1[jx][jy]; // X 1st derivative
    coef5 = mesh->G3[jx][jy]; // Z 1st derivative
  }

  if(d != (Field2D*) NULL) {
    // Multiply Delp2 component by a factor
    coef1 *= (*d)[jx][jy];
    coef2 *= (*d)[jx][jy];
    coef3 *= (*d)[jx][jy];
    coef4 *= (*d)[jx][jy];
    coef5 *= (*d)[jx][jy];
  }

  if(nonuniform) {
    // non-uniform mesh correction
    if((jx != 0) && (jx != (mesh->ngx-1))) {
      //coef4 += mesh->g11[jx][jy]*0.25*( (1.0/dx[jx+1][jy]) - (1.0/dx[jx-1][jy]) )/dx[jx][jy]; // SHOULD BE THIS (?)
      coef4 -= 0.5*((mesh->dx[jx+1][jy] - mesh->dx[jx-1][jy])/SQ(mesh->dx[jx][jy]))*coef1; // BOUT-06 term
    }
  }

  if(ccoef != NULL) {
    // A first order derivative term
    
    if((jx > 0) && (jx < (mesh->ngx-1)))
      coef4 += mesh->g11[jx][jy] * ((*ccoef)[jx+1][jy] - (*ccoef)[jx-1][jy]) / (2.*mesh->dx[jx][jy]*((*ccoef)[jx][jy]));
  }
  
  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // d2dz2 term
    coef2 += mesh->g11[jx][jy] * mesh->IntShiftTorsion[jx][jy] * mesh->IntShiftTorsion[jx][jy];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }
  
  coef1 /= SQ(mesh->dx[jx][jy]);
  coef3 /= 2.*mesh->dx[jx][jy];
  coef4 /= 2.*mesh->dx[jx][jy];

  a = dcomplex(coef1 - coef4,-kwave*coef3);
  b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,kwave*coef5);
  c = dcomplex(coef1 + coef4,kwave*coef3);
}

/// Sets the coefficients for parallel tridiagonal matrix inversion
/*!
 * Uses the laplace_tridag_coefs routine above to fill a matrix [kz][ix] of coefficients
 */
void Laplacian::tridagMatrix(dcomplex **avec, dcomplex **bvec, dcomplex **cvec,
                             dcomplex **bk, int jy, int flags, 
                             const Field2D *a, const Field2D *ccoef, 
                             const Field2D *d) {
  for(int kz = 0; kz <= maxmode; kz++) {
    BoutReal kwave=kz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
    
    
    tridagMatrix(avec[kz], bvec[kz], cvec[kz],
                 bk[kz], 
                 jy, 
                 kz == 0, kwave, 
                 flags, 
                 a, ccoef, d);
  }
}

void Laplacian::tridagMatrix(dcomplex *avec, dcomplex *bvec, dcomplex *cvec,
                             dcomplex *bk, int jy, bool dc, BoutReal kwave, 
                             int flags, 
                             const Field2D *a, const Field2D *ccoef, 
                             const Field2D *d,
                             bool includeguards) {
  int xs = 0;
  int xe = mesh->ngx-1;
  if(!includeguards) {
    if(!mesh->firstX()) 
      xs = mesh->xstart; // Inner edge is a guard cell
    if(!mesh->lastX())
      xe = mesh->xend; // Outer edge is a guard cell
  }
  int ncx = xe - xs;
  
  int inbndry = 2, outbndry=2;
  
  if(flags & INVERT_BNDRY_ONE) {
    inbndry = outbndry = 1;
  }
  if(flags & INVERT_BNDRY_IN_ONE)
    inbndry = 1;
  if(flags & INVERT_BNDRY_OUT_ONE)
    outbndry = 1;
    
  // Entire domain. Change to boundaries later
  
  for(int ix=0;ix<=ncx;ix++) {
    
    tridagCoefs(xs+ix, jy, kwave, avec[ix], bvec[ix], cvec[ix], ccoef, d);
      
    if(a != (Field2D*) NULL)
      bvec[ix] += (*a)[xs+ix][jy];
  }

  if(!mesh->periodicX) {
    // Boundary conditions
    
    if(mesh->firstX()) {
      // INNER BOUNDARY ON THIS PROCESSOR
      
      if(!(flags & (INVERT_IN_RHS | INVERT_IN_SET))) {
        for(int ix=0;ix<inbndry;ix++)
          bk[ix] = 0.;
      }
      
      if(dc) {
        // DC i.e. kz = 0
	  
        if(flags & INVERT_DC_IN_GRAD) {
          // Zero gradient at inner boundary
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] =  0.;
            bvec[ix] =  1.;
            cvec[ix] = -1.;
          }
        }else if(flags & INVERT_DC_IN_GRADPAR) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  1.0/sqrt(mesh->g_22(ix,jy));
            cvec[ix] = -1.0/sqrt(mesh->g_22(ix+1,jy));
          }
        }else if(flags & INVERT_DC_IN_GRADPARINV) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  sqrt(mesh->g_22(ix,jy));
            cvec[ix] = -sqrt(mesh->g_22(ix+1,jy));
          }
        }else if (flags & INVERT_DC_IN_LAP) {
          // Decaying boundary conditions
          BoutReal k = 0.0;
          if(a != (Field2D*) NULL) {
            BoutReal ksq = -((*a)(inbndry, jy));
            if(ksq < 0.0)
              throw BoutException("ksq must be positive");
            k = sqrt(ksq);
          }
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] =  0.;
            bvec[ix] =  1.;
            cvec[ix] = -exp(-k*mesh->dx(ix,jy)/sqrt(mesh->g11(ix,jy)));
          }
            
        }else {
          // Zero value at inner boundary or INVERT_IN_SET
          for (int ix=0;ix<inbndry;ix++){
            avec[ix] = 0.;
            bvec[ix] = 0.5;
            cvec[ix] = 0.5;
          }
        }
	  
      }else {
        // AC
	
        if(flags & INVERT_AC_IN_GRAD) {
          // Zero gradient at inner boundary
          for (int ix=0;ix<inbndry;ix++){
            avec[ix]=dcomplex(0.,0.);
            bvec[ix]=dcomplex(1.,0.);
            cvec[ix]=dcomplex(-1.,0.);
          }
        }else if(flags & INVERT_AC_IN_LAP) {
          // Use decaying zero-Laplacian solution in the boundary
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] = 0.0;
            bvec[ix] = 1.0;
            cvec[ix] = -exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*mesh->dx[ix][jy]);
          }
        }else {
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
      // OUTER BOUNDARY
      
      if(!(flags & (INVERT_OUT_RHS | INVERT_OUT_SET))) {
        for (int ix=0;ix<outbndry;ix++) {
          bk[ncx-ix] = 0.;
        }
      }

      if(dc) {
        // DC
	
        if(flags & INVERT_DC_OUT_GRAD) {
          // Zero gradient at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            cvec[ncx-ix]=dcomplex(0.,0.);
            bvec[ncx-ix]=dcomplex(1.,0.);
            avec[ncx-ix]=dcomplex(-1.,0.);
          }
        }else {
          // Zero value at outer boundary or INVERT_OUT_SET
          for (int ix=0;ix<outbndry;ix++){
            cvec[ncx-ix]=dcomplex(0.,0.);
            bvec[ncx-ix]=dcomplex(0.5,0.);
            avec[ncx-ix]=dcomplex(0.5,0.);
          }
        }
      }else {
        // AC
	
        if(flags & INVERT_AC_OUT_GRAD) {
          // Zero gradient at outer boundary
          for (int ix=0;ix<outbndry;ix++){
            cvec[ncx-ix]=dcomplex(0.,0.);
            bvec[ncx-ix]=dcomplex(1.,0.);
            avec[ncx-ix]=dcomplex(-1.,0.);
          }
        }else if(flags & INVERT_AC_OUT_LAP) {
          // Use decaying zero-Laplacian solution in the boundary
          for (int ix=0;ix<outbndry;ix++) {
            avec[ncx-ix] = -exp(-1.0*sqrt(mesh->g33[xe-ix][jy]/mesh->g11[xe-ix][jy])*kwave*mesh->dx[xe-ix][jy]);;
            bvec[ncx-ix] = 1.0;
            cvec[ncx-ix] = 0.0;
          }
        }else {
          // Zero value at outer boundary or LAPLACE_OUT_SET
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
                          const Field2D *ccoef, const Field2D *d) {
  Laplacian::defaultInstance()->tridagCoefs(jx,jy, jz, a, b, c, ccoef, d);
}

int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  
  Laplacian *lap = Laplacian::defaultInstance();
  
  if(a != NULL) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);
  
  if(c != NULL) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);
  
  if(d != NULL) {
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
  
  if(a != NULL) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);
  
  if(c != NULL) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);
  
  if(d != NULL) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);
  
  lap->setFlags(flags);
  
  x.allocate(); // Make sure x is allocated

  x = lap->solve(b, x);
  
  x.setLocation(b.getLocation());

}
const Field3D invert_laplace(const Field3D &b, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  
  Timer timer("invert"); ///< Start timer 
  
  Laplacian *lap = Laplacian::defaultInstance();
  
  if(a != NULL) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);
  
  if(c != NULL) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);
  
  if(d != NULL) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);
  
  lap->setFlags(flags);
  
  Field3D x = lap->solve(b);
  
  x.setLocation(b.getLocation());

  return x;
}

