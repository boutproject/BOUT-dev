/**************************************************************************
 * Basic derivative methods
 *
 * 
 * Four kinds of differencing methods:
 * 
 * 1. First derivative DD*
 *    Central differencing e.g. Div(f)
 *
 * 2. Second derivatives D2D*2
 *    Central differencing e.g. Delp2(f)
 *
 * 3. Upwinding VDD*
 *    Terms like v*Grad(f)
 *
 * 4. Flux methods FDD* (e.g. flux conserving, limiting)
 *    Div(v*f)
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
 **************************************************************************/

#include <globals.hxx>
#include <derivs.hxx>
#include <stencils.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <interpolation.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>

#include <cmath>
#include <string.h>
#include <stdlib.h>

#include <output.hxx>

//#undef _OPENMP

#ifdef _OPENMP
#include <omp.h>
#endif

/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D DDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->indexDDX(f,outloc, method) / mesh->coordinates()->dx;
}

const Field3D DDX(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return DDX(f, outloc, method);
}

const Field3D DDX(const Field3D &f, DIFF_METHOD method) {
  return DDX(f, CELL_DEFAULT, method);
}

const Field2D DDX(const Field2D &f) {
  return mesh->indexDDX(f) / mesh->coordinates()->dx;
}

////////////// Y DERIVATIVE /////////////////

const Field3D DDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->indexDDY(f,outloc, method) / mesh->coordinates()->dy;
}

const Field3D DDY(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return DDY(f, outloc, method);
}

const Field3D DDY(const Field3D &f, DIFF_METHOD method) {
  return DDY(f, CELL_DEFAULT, method);
}

const Field2D DDY(const Field2D &f) {
  return mesh->indexDDY(f) / mesh->coordinates()->dy;
}

/*
const Field3D DDY_MUSCL(const Field3D &F, const Field3D &u, const Field2D &Vmax) {
  Field3D result;
  result.allocate(); // Make sure data allocated
  
  bindex bx;
  start_index(&bx, RGN_NOBNDRY);

  stencil fs, us;
  do {
    for(bx.jz=0;bx.jz<mesh->ngz-1;bx.jz++) {
      F.setYStencil(fs, bx);
      u.setYStencil(us, bx);
      
      result(bx.jx,bx.jy,bx.jz) = DDX_KT(fs, us, Vmax(bx.jx,bx.jy)) / mesh->dy(bx.jx, bx.jy);
    }
  }while(next_index2(&bx));
  
  return result;
}
*/

////////////// Z DERIVATIVE /////////////////

const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, bool inc_xbndry) {
  return mesh->indexDDZ(f,outloc, method, inc_xbndry) / mesh->coordinates()->dz;
}

const Field3D DDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc, bool inc_xbndry) {
  return DDZ(f, outloc, method, inc_xbndry);
}

const Field3D DDZ(const Field3D &f, DIFF_METHOD method, bool inc_xbndry) {
  return DDZ(f, CELL_DEFAULT, method, inc_xbndry);
}

const Field3D DDZ(const Field3D &f, bool inc_xbndry) {
  return DDZ(f, CELL_DEFAULT, DIFF_DEFAULT, inc_xbndry);
}

const Field2D DDZ(const Field2D &f) {
  Field2D result;
  result = 0.0;
  return result;
}

const Vector3D DDZ(const Vector3D &v, CELL_LOC outloc, DIFF_METHOD method) {
  Vector3D result;

  result.covariant = v.covariant;

  result.x = DDZ(v.x, outloc, method);
  result.y = DDZ(v.y, outloc, method);
  result.z = DDZ(v.z, outloc, method);

  return result;
}

const Vector3D DDZ(const Vector3D &v, DIFF_METHOD method, CELL_LOC outloc) {
  return DDZ(v, outloc, method);
}

const Vector2D DDZ(const Vector2D &v) {
  Vector2D result;

  result.covariant = v.covariant;

  result.x = 0.;
  result.y = 0.;
  result.z = 0.;

  return result;
}

/*******************************************************************************
 * 2nd derivative
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D D2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  deriv_func func = fD2DX2; // Set to default function
  inner_boundary_deriv_func func_in = fD2DX2_in;
  outer_boundary_deriv_func func_out = fD2DX2_out;
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) ||
       ((inloc == CELL_XLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in X. Centre -> Xlow, or Xlow -> Centre
      
      func = sfD2DX2; // Set default
      func_in = sfD2DX2_in;
      func_out = sfD2DX2_out;
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_XLOW) {
	// Shifting
	
	func = sfD2DX2; // Set default
	func_in = sfD2DX2_in;
	func_out = sfD2DX2_out;
	table = SecondStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return D2DX2(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    func_in = lookupInnerBoundaryFunc(table, method);
    func_out = lookupOuterBoundaryFunc(table, method);
    if(func == NULL)
      bout_error("Cannot use FFT for X derivatives");
  }
  
  Field2D dd = mesh->dx*mesh->dx;
  result = applyXdiff(f, func, func_in, func_out, dd);
  result.setLocation(diffloc);
  
  if(non_uniform) {
    // Correction for non-uniform mesh
    result += mesh->d1_dx*applyXdiff(f, fDDX, fDDX_in, fDDX_out, mesh->dx);
  }
  
  result = interp_to(result, outloc);

  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    mesh->IncIntShear = false; // So DDX doesn't try to include I again
    // Add I^2 d^2/dz^2 term
    result += mesh->IntShiftTorsion^2 * D2DZ2(f, outloc);
    // Mixed derivative
    result += 2.*mesh->IntShiftTorsion * D2DXDZ(f);
    // DDZ term
    result += DDX(mesh->IntShiftTorsion) * DDZ(f, outloc);
    mesh->IncIntShear = true;
  }
  
  return result;
}

const Field3D D2DX2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return D2DX2(f, outloc, method);
}

const Field2D D2DX2(const Field2D &f) {
  Field2D result;

  Field2D dd = mesh->dx*mesh->dx;
  result = applyXdiff(f, fD2DX2, fD2DX2_in, fD2DX2_out, dd);
  
  if(non_uniform) {
    // Correction for non-uniform mesh
    result += mesh->d1_dx * applyXdiff(f, fDDX, fDDX_in, fDDX_out, mesh->dx);
  }
  
  return(result);
}

////////////// Y DERIVATIVE /////////////////

const Field3D D2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  deriv_func func = fD2DY2; // Set to default function
  inner_boundary_deriv_func func_in = fD2DY2_in;
  outer_boundary_deriv_func func_out = fD2DY2_out;
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_YLOW)) ||
       ((inloc == CELL_YLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Y. Centre -> Ylow, or Ylow -> Centre
      
      func = sfD2DY2; // Set default
      func_in = sfD2DY2_in;
      func_out = sfD2DY2_out;
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_YLOW) {
	// Shifting
	
	func = sfD2DY2; // Set default
	func_in = sfD2DY2_in;
	func_out = sfD2DY2_out;
	table = SecondStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return D2DY2(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    func_in = lookupInnerBoundaryFunc(table, method);
    func_out = lookupOuterBoundaryFunc(table, method);
    if(func == NULL)
      bout_error("Cannot use FFT for Y derivatives");
  }
  
  Field2D dd = mesh->dy*mesh->dy;
  result = applyYdiff(f, func, func_in, func_out, dd);
  result.setLocation(diffloc);

  if(non_uniform) {
    // Correction for non-uniform mesh
    result += mesh->d1_dy * applyYdiff(f, fDDY, fDDY_in, fDDY_out, mesh->dy);
  }

  return interp_to(result, outloc);
}

const Field3D D2DY2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return D2DY2(f, outloc, method);
}

const Field2D D2DY2(const Field2D &f) {
  return applyYdiff(f, fD2DY2, fD2DY2_in, fD2DY2_out, mesh->dy*mesh->dy);
}

////////////// Z DERIVATIVE /////////////////

const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  deriv_func func = fD2DZ2; // Set to default function
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;

  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_ZLOW)) ||
       ((inloc == CELL_ZLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Z. Centre -> Zlow, or Zlow -> Centre
      
      func = sfD2DZ2; // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_ZLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_ZLOW) {
	// Shifting
	
	func = sfD2DZ2; // Set default
	table = SecondStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return D2DZ2(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if(func == NULL) {
    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?
    if(mesh->StaggerGrids) {
      if((inloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
	// Shifting down - multiply by exp(-0.5*i*k*dz)
	shift = -1.;
      }else if((inloc == CELL_ZLOW) && (diffloc == CELL_CENTRE)) {
	// Shifting up
	shift = 1.;
      }
    }
    
    result.allocate(); // Make sure data allocated

    int ncz = mesh->ngz-1;
    
#ifndef _OPENMP
    static dcomplex *cv = (dcomplex*) NULL;
#else
    static dcomplex *globalcv;
    static int nthreads = 0;
#endif
    
    #pragma omp parallel
    {
#ifndef _OPENMP
      // Serial, so can have a single static array
      if(cv == (dcomplex*) NULL)
        cv = new dcomplex[ncz/2 + 1];
#else
      // Parallel, so allocate a separate array for each thread
      
      int th_id = omp_get_thread_num(); // thread ID
      int n_th = omp_get_num_threads();
      if(th_id == 0) {
        if(nthreads < n_th) {
          // Allocate memory in thread zero
          if(nthreads > 0)
            delete[] globalcv;
          globalcv = new dcomplex[n_th*(ncz/2 + 1)];
          nthreads = n_th;
        }
      }
      // Wait for memory to be allocated
      #pragma omp barrier
      
      dcomplex *cv = globalcv + th_id*(ncz/2 + 1); // Separate array for each thread
#endif
      int xs = mesh->xstart;
      int xe = mesh->xend;
      int ys = mesh->ystart;
      int ye = mesh->yend;
      if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX)
	xs = 0;
      if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX)
	xe = mesh->ngx-1;
      if (mesh->freeboundary_ydown)
	ys = 0;
      if (mesh->freeboundary_ydown)
	ye = mesh->ngy-1;
      #pragma omp for
      for(int jx=xs;jx<=xe;jx++) {
        for(int jy=ys;jy<=ye;jy++) {
          
          rfft(f[jx][jy], ncz, cv); // Forward FFT
	
          for(int jz=0;jz<=ncz/2;jz++) {
            BoutReal kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
            
            BoutReal flt;
            if (jz>0.4*ncz) flt=1e-10; else flt=1.0;

            cv[jz] *= -SQ(kwave) * flt;
            if(mesh->StaggerGrids)
              cv[jz] *= exp(Im * (shift * kwave * mesh->dz));
          }

          irfft(cv, ncz, result[jx][jy]); // Reverse FFT
	
          result(jx,jy,ncz) = result(jx,jy,0);
        }
      }
    } // End of parallel section

#ifdef CHECK
    // Mark boundaries as invalid
    if (mesh->freeboundary_xin) result.bndry_xin = true;
    else result.bndry_xin = false;
    if (mesh->freeboundary_xout) result.bndry_xout = true;
    else result.bndry_xout = false;
    if (mesh->freeboundary_yup) result.bndry_yup = true;
    else result.bndry_yup = false;
    if (mesh->freeboundary_ydown) result.bndry_ydown = true;
    else result.bndry_ydown = false;
#endif

  }
  else {
    // All other (non-FFT) functions
    result = applyZdiff(f, func, SQ(mesh->dz));
  }

  result.setLocation(diffloc);

  return interp_to(result, outloc);
}

const Field3D D2DZ2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return D2DZ2(f, outloc, method);
}


const Field2D D2DZ2(const Field2D &f) {
  Field2D result;
  result = 0.0;
  return result;
}

/*******************************************************************************
 * Fourth derivatives
 *******************************************************************************/

BoutReal D4DX4_C2(stencil &f) {
  return (f.pp - 4.*f.p + 6.*f.c - 4.*f.m + f.mm);
}

boundary_derivs_pair D4D4_F2(forward_stencil &f) {
  boundary_derivs_pair result;
  result.inner = 2.*f.m-9.*f.c+16.*f.p-14.*f.p2+6.*f.p3-f.p4;
  result.outer = 3.*f.m-14.*f.c+26.*f.p-24.*f.p2+11.*f.p3-2.*f.p4;
}

boundary_derivs_pair D4D4_B2(backward_stencil &f) {
  boundary_derivs_pair result;
  result.inner = 2.*f.p-9.*f.c+16.*f.m-14.*f.m2+6.*f.m3-f.m4;
  result.outer = 3.*f.p-14.*f.c+26.*f.m-24.*f.m2+11.*f.m3-2.*f.m4;
}

const Field3D D4DX4(const Field3D &f) {
  return applyXdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2, SQ(SQ(mesh->dx)));
}

const Field2D D4DX4(const Field2D &f) {
  return applyXdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2, SQ(SQ(mesh->dx)));
}

const Field3D D4DY4(const Field3D &f) {
  return applyYdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2, SQ(SQ(mesh->dy)));
}

const Field2D D4DY4(const Field2D &f) {
  return applyYdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2, SQ(SQ(mesh->dy)));
}

const Field3D D4DZ4(const Field3D &f) {
  return applyZdiff(f, D4DX4_C2, SQ(SQ(mesh->dz)));
}

const Field2D D4DZ4(const Field2D &f) {
  return Field2D(0.0);
}

/*******************************************************************************
 * Mixed derivatives
 *******************************************************************************/

const Field2D D2DXDY(const Field2D &f) {
  return Field2D(0.0);
  /*
    Note: Missing corners, so following will break
  Field2D result;
  result.allocate();
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      result(i,j) = 0.25*( +(f(i+1,j+1) - f(i-1,j+1))/(mesh->dx(i,j+1))
                           -(f(i+1,j-1) - f(i-1,j-1))/(mesh->dx(i,j-1)) )
        / mesh->dy(i,j);
    }
  return result;
  */
}

const Field3D D2DXDY(const Field3D &f) {
  return Field3D(0.0);
  /*
    Note: Missing corners, so following will break
    
  Field3D result;
  result.allocate();
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) 
      for(int k=0;k<mesh->ngz;k++) {
        result(i,j,k) = 0.25*( +(f(i+1,j+1,k) - f(i-1,j+1,k))/(mesh->dx(i,j+1))
                               -(f(i+1,j-1,k) - f(i-1,j-1,k))/(mesh->dx(i,j-1)) )
          / mesh->dy(i,j);
      }
  return result;
  */
}

const Field2D D2DXDZ(const Field2D &f) {
  return Field2D(0.0);
}

/// X-Z mixed derivative
const Field3D D2DXDZ(const Field3D &f) {
  Field3D result;
  
  // Take derivative in Z, including in X boundaries. Then take derivative in X
  // Maybe should average results of DDX(DDZ) and DDZ(DDX)?
  result = DDX(DDZ(f, true));
  
  return result;
}

const Field2D D2DYDZ(const Field2D &f) {
  return Field2D(0.0);
}

const Field3D D2DYDZ(const Field3D &f) {
  Field3D result;
  result.allocate();
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) 
      for(int k=0;k<mesh->ngz-1;k++) {
        int kp = (k+1) % (mesh->ngz-1);
        int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
        result(i,j,k) = 0.25*( +(f(i,j+1,kp) - f(i,j-1,kp))/(mesh->dy(i,j+1))
                               -(f(i,j+1,km) - f(i,j-1,km))/(mesh->dy(i,j-1)) )
          / mesh->dz;
      }
  return result;
}

/*******************************************************************************
 * Advection schemes
 * 
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/// Special case where both arguments are 2D. Output location ignored for now
const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  upwind_func func = fVDDX;

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(UpwindTable, method);
  }

  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setXStencil(fs, bx);
    v.setXStencil(vs, bx);
    
    d[bx.jx][bx.jy] = func(vs, fs) / mesh->dx(bx.jx, bx.jy);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  return result;
}

const Field2D VDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method) {
  return VDDX(v, f, CELL_DEFAULT, method);
}

/// General version for 2 or 3-D objects
const Field3D VDDX(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  upwind_func func = fVDDX;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_XLOW) {
      // V staggered w.r.t. variable
      func = sfVDDX;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_XLOW)) {
      // Shifted
      func = sfVDDX;
      table = UpwindStagTable;
      diffloc = CELL_XLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  /// Clone inputs (for shifting)
  Field *vp = v.clone();
  Field *fp = f.clone();
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    // Shift in Z using FFT if needed
    vp->shiftToReal(true);
    fp->shiftToReal(true);
  }
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  bindex bx;
  
  start_index(&bx);
#ifdef _OPENMP
  // Parallel version
  
  bindex bxstart = bx; // Copy to avoid race condition on first index
  bool workToDoGlobal; // Shared loop control
  #pragma omp parallel
  {
    bindex bxlocal; // Index for each thread
    stencil vval, fval;
    bool workToDo;  // Does this thread have work to do?
    #pragma omp single
    {
      // First index done by a single thread
      for(bxstart.jz=0;bxstart.jz<mesh->ngz-1;bxstart.jz++) {
        vp->setXStencil(vval, bxstart, diffloc);
        fp->setXStencil(fval, bxstart); // Location is always the same as input
    
        d[bxstart.jx][bxstart.jy][bxstart.jz] = func(vval, fval) / mesh->dx(bxstart.jx, bxstart.jy);
      }
    }
    
    do {
      #pragma omp critical
      {
        // Get the next index
        workToDo = next_index2(&bx);
        bxlocal = bx; // Make a local copy
        workToDoGlobal = workToDo;
      }
      if(workToDo) { // Here workToDo could be different to workToDoGlobal
        for(bxlocal.jz=0;bxlocal.jz<mesh->ngz-1;bxlocal.jz++) {
          vp->setXStencil(vval, bxlocal, diffloc);
          fp->setXStencil(fval, bxlocal); // Location is always the same as input
    
          d[bxlocal.jx][bxlocal.jy][bxlocal.jz] = func(vval, fval) / mesh->dx(bxlocal.jx, bxlocal.jy);
        }
      }
    }while(workToDoGlobal);
  }
#else
  // Serial version
  stencil vval, fval;
  do {
    vp->setXStencil(vval, bx, diffloc);
    fp->setXStencil(fval, bx); // Location is always the same as input
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval) / mesh->dx(bx.jx, bx.jy);
  }while(next_index3(&bx));
#endif
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
    result = result.shiftZ(false); // Shift back
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  // Delete clones
  delete vp;
  delete fp;
  
  return interp_to(result, outloc);
}

const Field3D VDDX(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc) {
  return VDDX(v, f, outloc, method);
}

////////////// Y DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  upwind_func func = fVDDY;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_YLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDY(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  bindex bx;
  stencil fval, vval;
  
  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  start_index(&bx);
  do {
    f.setYStencil(fval, bx);
    v.setYStencil(vval, bx, diffloc);
    d[bx.jx][bx.jy] = func(vval,fval)/mesh->dy(bx.jx, bx.jy);
  }while(next_index2(&bx));

  result.setLocation(inloc);
  
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

const Field2D VDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method) {
  return VDDY(v, f, CELL_DEFAULT, method);
}

// general case
const Field3D VDDY(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  upwind_func func = fVDDY;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_YLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDY(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }
  
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }
  
  bindex bx;
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  start_index(&bx);
#ifdef _OPENMP
  // Parallel version
  
  bindex bxstart = bx; // Copy to avoid race condition on first index
  bool workToDoGlobal; // Shared loop control

  #pragma omp parallel
  {
    bindex bxlocal; // Index for each thread
    stencil vval, fval;
    bool workToDo;  // Does this thread have work to do?
    #pragma omp single
    {
      // First index done by a single thread
      for(bxstart.jz=0;bxstart.jz<mesh->ngz-1;bxstart.jz++) {
        v.setYStencil(vval, bxstart, diffloc);
        f.setYStencil(fval, bxstart);
    
        d[bxstart.jx][bxstart.jy][bxstart.jz] = func(vval, fval)/mesh->dy(bxstart.jx, bxstart.jy);
      }
    }
    
    do {
      #pragma omp critical
      {
        // Get the next index
        workToDo = next_index2(&bx);
        bxlocal = bx; // Make a local copy
        workToDoGlobal = workToDo;
      }
      if(workToDo) { // Here workToDo could be different to workToDoGlobal
        for(bxlocal.jz=0;bxlocal.jz<mesh->ngz-1;bxlocal.jz++) {
          v.setYStencil(vval, bxlocal, diffloc);
          f.setYStencil(fval, bxlocal);
    
          d[bxlocal.jx][bxlocal.jy][bxlocal.jz] = func(vval, fval)/mesh->dy(bxlocal.jx, bxlocal.jy);
        }
      }
    }while(workToDoGlobal);
  }
#else
  // Serial version
  stencil vval, fval;
  do {
    v.setYStencil(vval, bx, diffloc);
    f.setYStencil(fval, bx);
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval)/mesh->dy(bx.jx, bx.jy);
  }while(next_index3(&bx));
#endif
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

const Field3D VDDY(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc) {
  return VDDY(v, f, outloc, method);
}

////////////// Z DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDZ(const Field2D &v, const Field2D &f) {
  Field2D result;
  result = 0.0;
  return result;
}

// Note that this is zero because no compression is included
const Field2D VDDZ(const Field3D &v, const Field2D &f) {
  Field2D result;
  result = 0.0;
  return result;
}

// general case
const Field3D VDDZ(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  upwind_func func = fVDDZ;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_ZLOW) {
      // V staggered w.r.t. variable
      func = sfVDDZ;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_ZLOW)) {
      // Shifted
      func = sfVDDZ;
      table = UpwindStagTable;
      diffloc = CELL_ZLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDY(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  bindex bx;
  
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();
  
  start_index(&bx);
#ifdef _OPENMP
  // Parallel version

  bindex bxstart = bx; // Copy to avoid race condition on first index
  bool workToDoGlobal; // Shared loop control
  
  #pragma omp parallel
  {
    bindex bxlocal; // Index for each thread
    stencil vval, fval;
    bool workToDo;  // Does this thread have work to do?
    #pragma omp single
    {
      // First index done by a single thread
      for(bxstart.jz=0;bxstart.jz<mesh->ngz-1;bxstart.jz++) {
        v.setZStencil(vval, bxstart, diffloc);
        f.setZStencil(fval, bxstart);
    
        d[bxstart.jx][bxstart.jy][bxstart.jz] = func(vval, fval)/mesh->dz;
      }
    }
    
    do {
      #pragma omp critical
      {
        // Get the next index
        workToDo = next_index2(&bx);
        bxlocal = bx; // Make a local copy
        workToDoGlobal = workToDo;
      }
      if(workToDo) { // Here workToDo could be different to workToDoGlobal
        for(bxlocal.jz=0;bxlocal.jz<mesh->ngz-1;bxlocal.jz++) {
          v.setZStencil(vval, bxlocal, diffloc);
          f.setZStencil(fval, bxlocal);
    
          d[bxlocal.jx][bxlocal.jy][bxlocal.jz] = func(vval, fval)/mesh->dz;
        }
      }
    }while(workToDoGlobal);
  }
#else 
  stencil vval, fval;
  do {
    v.setZStencil(vval, bx, diffloc);
    f.setZStencil(fval, bx);
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval)/mesh->dz;
  }while(next_index3(&bx));
#endif

  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

const Field3D VDDZ(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc) {
  return VDDZ(v, f, outloc, method);
}

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/

const Field2D FDDX(const Field2D &v, const Field2D &f) {
  return FDDX(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

const Field2D FDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return FDDX(v, f, method, outloc);
}

const Field2D FDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method, CELL_LOC outloc) {
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return VDDX(v, f) + f * DDX(v);
  }
  
  upwind_func func = fFDDX;
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(FluxTable, method);
  }
  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setXStencil(fs, bx);
    v.setXStencil(vs, bx);
    
    d[bx.jx][bx.jy] = func(vs, fs) / mesh->dx(bx.jx, bx.jy);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif
  
  return result;
}

const Field3D FDDX(const Field3D &v, const Field3D &f) {
  return FDDX(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

const Field3D FDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return FDDX(v, f, method, outloc);
}

const Field3D FDDX(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return VDDX(v, f, outloc) + DDX(v, outloc) * f;
  }
  
  upwind_func func = fFDDX;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_XLOW) {
      // V staggered w.r.t. variable
      func = sfFDDX;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_XLOW)) {
      // Shifted
      func = sfFDDX;
      table = FluxStagTable;
      diffloc = CELL_XLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  /// Clone inputs (for shifting)
  Field3D *vp = v.clone();
  Field3D *fp = f.clone();
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    // Shift in Z using FFT if needed
    vp->shiftToReal(true);
    fp->shiftToReal(true);
  }
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    vp->setXStencil(vval, bx, diffloc);
    fp->setXStencil(fval, bx); // Location is always the same as input
    
    result(bx.jx,bx.jy,bx.jz) = func(vval, fval) / mesh->dx(bx.jx,bx.jy);
  }while(next_index3(&bx));
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
    result = result.shiftZ(false); // Shift back
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  // Delete clones
  delete vp;
  delete fp;
  
  return interp_to(result, outloc);
}

/////////////////////////////////////////////////////////////////////////

const Field2D FDDY(const Field2D &v, const Field2D &f) {
  return FDDY(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

const Field2D FDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return FDDY(v, f, method, outloc);
}

const Field2D FDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method, CELL_LOC outloc) {
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return VDDY(v, f) + f * DDY(v);
  }
  
  upwind_func func = fFDDY;
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(FluxTable, method);
  }
  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setYStencil(fs, bx);
    v.setYStencil(vs, bx);
    
    d[bx.jx][bx.jy] = func(vs, fs) / mesh->dy(bx.jx, bx.jy);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif
  
  return result;
}

const Field3D FDDY(const Field3D &v, const Field3D &f) {
  return FDDY(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

const Field3D FDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return FDDY(v, f, method, outloc);
}

const Field3D FDDY(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  //output.write("fddy: %d %d %d : %u %u \n", v.getLocation(), f.getLocation(), method, fFDDY, sfFDDY);
  
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return VDDY(v, f, outloc) + DDY(v, outloc) * f;
  }
  upwind_func func = fFDDY;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value    
    if(vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfFDDY;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfFDDY;
      table = FluxStagTable;
      diffloc = CELL_YLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.
      
      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }
  
  if(func == NULL) {
    // To catch when no function
    return VDDY(v, f, outloc) + DDY(v, outloc) * f;
  }

  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    v.setYStencil(vval, bx, diffloc);
    f.setYStencil(fval, bx); // Location is always the same as input
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval) / mesh->dy(bx.jx, bx.jy);

  }while(next_index3(&bx));

  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

/////////////////////////////////////////////////////////////////////////

const Field2D FDDZ(const Field2D &v, const Field2D &f) {
  return FDDZ(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

const Field2D FDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return FDDZ(v, f, method, outloc);
}

const Field2D FDDZ(const Field2D &v, const Field2D &f, DIFF_METHOD method, CELL_LOC outloc) {
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDZ == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return VDDZ(v, f) + f * DDZ(v);
  }
  
  upwind_func func = fFDDZ;
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(FluxTable, method);
  }
  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setZStencil(fs, bx);
    v.setZStencil(vs, bx);
    
    d[bx.jx][bx.jy] = func(vs, fs) / mesh->dz;
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif
  
  return result;
}

const Field3D FDDZ(const Field3D &v, const Field3D &f) {
  return FDDZ(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

const Field3D FDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return FDDZ(v, f, method, outloc);
}

const Field3D FDDZ(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDZ == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return VDDZ(v, f, outloc) + DDZ(v, outloc) * f;
  }
  
  upwind_func func = fFDDZ;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_ZLOW) {
      // V staggered w.r.t. variable
      func = sfFDDZ;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_ZLOW)) {
      // Shifted
      func = sfFDDZ;
      table = FluxStagTable;
      diffloc = CELL_ZLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    v.setZStencil(vval, bx, diffloc);
    f.setZStencil(fval, bx); // Location is always the same as input
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval) / mesh->dz;
  }while(next_index3(&bx));
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  return interp_to(result, outloc);
}
