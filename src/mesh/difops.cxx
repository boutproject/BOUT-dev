/**************************************************************************
 * Various differential operators defined on BOUT grid
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
#include <bout.hxx>
#include <difops.hxx>
#include <vecops.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <msg_stack.hxx>
#include <bout/assert.hxx>

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

#include <interpolation.hxx>

#include <math.h>
#include <stdlib.h>

/*******************************************************************************
 * Grad_par
 * The parallel derivative along unperturbed B-field
 *******************************************************************************/

const Field2D Grad_par(const Field2D &var, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->coordinates()->Grad_par(var, outloc, method);
}

const Field2D Grad_par(const Field2D &var, DIFF_METHOD method, CELL_LOC outloc) {
  return mesh->coordinates()->Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->coordinates()->Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, DIFF_METHOD method, CELL_LOC outloc) {
  return mesh->coordinates()->Grad_par(var, outloc, method);
}

// Model dvar/dt = Grad_par(f) with a maximum velocity of Vmax
const Field3D Grad_par(const Field3D &f, const Field3D &var, const Field2D &Vmax) {
  int msg_pos = msg_stack.push("Grad_par( Field3D, Field3D, Field2D )");

  Field2D sg = sqrt(mesh->coordinates()->g_22);
  Field3D result = DDY_MUSCL(f, var, sg*Vmax)/sg;
  
  msg_stack.pop(msg_pos);
  
  return result;
}

const Field3D Grad_par(const Field3D &f, const Field3D &var, BoutReal Vmax) {
  Field2D V = Vmax;
  return Grad_par(f, var, V);
}

/*******************************************************************************
 * Grad_parP
 *
 * Derivative along perturbed field-line
 *
 * b0 dot Grad  -  (1/B)b0 x Grad(apar) dot Grad
 *
 * Combines the parallel and perpendicular calculation to include
 * grid-points at the corners.
 *******************************************************************************/

const Field3D Grad_parP(const Field3D &apar, const Field3D &f) {
  Field3D result;
  result.allocate();
  
  int ncz = mesh->ngz-1;

  Coordinates *metric = mesh->coordinates();
  
  Field3D gys;
  gys.allocate();

  // Need Y derivative everywhere
  for(int x=1;x<=mesh->ngx-2;x++)
    for(int y=1;y<=mesh->ngy-2;y++)
      for(int z=0;z<ncz;z++) {
        gys(x, y, z) = (f(x, y+1, z) - f(x, y-1, z))/(0.5*metric->dy(x, y+1) + metric->dy(x, y) + 0.5*metric->dy(x, y-1));
      }

  // Shift into orthogonal XZ local coordinates
  Field3D as = apar;
  Field3D fs = f;
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    as = apar.shiftZ(true);
    fs = f.shiftZ(true);
    gys = gys.shiftZ(true);
  }
  
  Field3D gx, bx, bz;
  gx.allocate();
  bx.allocate();
  bz.allocate();
  
  for(int x=1;x<=mesh->ngx-2;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      BoutReal by = 1./sqrt(metric->g_22(x, y));
      for(int z=0;z<ncz;z++) {
        int zm = (z - 1 + ncz) % ncz;
        int zp = (z + 1) % ncz;
        
        // bx = -DDZ(apar)
        bx(x, y, z) = (as(x, y, zm) - as(x, y, zp))/(2.*metric->dz);
        // bz = DDX(f)
        bz(x, y, z) = (as(x+1, y, z) - as(x-1, y, z))/(0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y));
        
        // Now calculate (bx*d/dx + by*d/dy + bz*d/dz) f
        
        // Length dl for predictor
        BoutReal dl = fabs(metric->dx(x, y)) / (fabs(bx(x, y, z)) + 1e-16);
        dl = BOUTMIN(dl, fabs(metric->dy(x, y)) / (fabs(by) + 1e-16));
        dl = BOUTMIN(dl, metric->dz / (fabs(bz(x, y, z)) + 1e-16));
        
        BoutReal fp, fm;
        
        // X differencing
        fp = fs(x+1, y, z)
          + (0.25*dl/metric->dz) * bz(x, y, z) * (fs(x+1, y, zm) - fs(x+1, y, zp))
          - 0.5*dl * by * gys(x+1, y, z);
        
        fm = fs(x-1, y, z)
          + (0.25*dl/metric->dz) * bz(x, y, z) * (fs(x-1, y, zm) - fs(x-1, y, zp))
          - 0.5*dl * by * gys(x-1, y, z);
        
        result(x, y, z) = bx(x, y, z) * (fp - fm) / (0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y));

        // Z differencing
        
        fp = fs(x, y, zp)
          + (0.25*dl/metric->dx(x, y)) * bx(x, y, z) * (fs(x-1, y, zp) - fs(x+1, y, zp))
          - 0.5*dl * by * gys(x, y, zp);
        
        fm = fs(x, y, zm)
          + (0.25*dl/metric->dx(x, y)) * bx(x, y, z) * (fs[x-1][y][zm] - fs(x+1, y, zm))
          - 0.5*dl * by * gys(x, y, zm);

        result(x, y, z) += bz(x, y, z) * (fp - fm) / (2.*metric->dz);
        
        // Y differencing. Need X derivative
        gx(x, y, z) = (fs(x+1, y, z) - fs(x-1, y, z))/(0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y));
        
      }
    }
  }

  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    // Shift into field-aligned coordinates
    gx = gx.shiftZ(false);
    bx = bx.shiftZ(false);
    bz = bz.shiftZ(false);
    result = result.shiftZ(false);
  }
  
  // Y differencing
  for(int x=1;x<=mesh->ngx-2;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      BoutReal by = 1./sqrt(metric->g_22(x, y));
      for(int z=0;z<ncz;z++) {
        int zm = (z - 1 + ncz) % ncz;
        int zp = (z + 1) % ncz;
        
        // Re-calculate dl
        BoutReal dl = fabs(metric->dx(x, y)) / (fabs(bx(x, y, z)) + 1e-16);
        dl = BOUTMIN(dl, fabs(metric->dy(x, y)) / (fabs(by) + 1e-16));
        dl = BOUTMIN(dl, metric->dz / (fabs(bz(x, y, z)) + 1e-16));
        
        BoutReal fp, fm;

        fp = f[x][y+1][z]
          - 0.5*dl * bx[x][y][z] * gx[x][y+1][z]
          + (0.25*dl/metric->dz)   * bz[x][y][z] * (fs[x][y+1][zm] - fs[x][y+1][zp]);
        
        fm = f[x][y-1][z]
          - 0.5*dl * bx[x][y][z] * gx[x][y-1][z]
          + (0.25*dl/metric->dz)   * bz[x][y][z] * (fs[x][y-1][zm] - fs[x][y-1][zp]);

        result[x][y][z] += by * (fp - fm) / (0.5*metric->dy[x][y-1] + metric->dy[x][y] + 0.5*metric->dy[x][y+1]);
      }
    }
  }
  
  return result;
}

/*******************************************************************************
 * Vpar_Grad_par
 * vparallel times the parallel derivative along unperturbed B-field
 *******************************************************************************/

const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f) {
  return mesh->coordinates()->Vpar_Grad_par(v, f);
}

const Field3D Vpar_Grad_par(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->coordinates()->Vpar_Grad_par(v, f, outloc, method);
}

const Field3D Vpar_Grad_par(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc) {
  return mesh->coordinates()->Vpar_Grad_par(v, f, outloc, method);
}

/*******************************************************************************
 * Div_par
 * parallel divergence operator B \partial_{||} (F/B)
 *******************************************************************************/

const Field2D Div_par(const Field2D &f) {
  return mesh->coordinates()->Div_par(f);
}

const Field3D Div_par(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->coordinates()->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return mesh->coordinates()->Div_par(f, outloc, method);
}

//////// Flux methods

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  Coordinates *metric = mesh->coordinates();
  return -metric->Bxy*FDDY(v, f/metric->Bxy, outloc, method)/sqrt(metric->g_22);
}

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return Div_par_flux(v,f, outloc, method);
}

//////// MUSCL schemes

const Field3D Div_par(const Field3D &f, const Field3D &var, const Field2D &Vmax) {
#ifdef CHECK
  int msg_pos = msg_stack.push("Div_par( Field3D, Field3D, Field2D )");
#endif
  Coordinates *metric = mesh->coordinates();
  
  Field3D result = metric->Bxy*Grad_par(f/metric->Bxy, var, Vmax/metric->Bxy);
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D Div_par(const Field3D &f, const Field3D &var, BoutReal Vmax) {
  Field2D V = Vmax;
  return Div_par(f, var, V);
}

/*******************************************************************************
 * Parallel derivatives converting between left and cell centred
 * NOTE: These are a quick hack to test if this works. The whole staggered grid
 *       thing needs to be thought through.
 *******************************************************************************/

const Field3D Grad_par_CtoL(const Field3D &var) {
  Field3D result;
  result.allocate();
  
  Coordinates *metric = mesh->coordinates();
  
  // NOTE: Need to calculate one more point than centred vars
  for(int jx=0; jx<mesh->ngx;jx++) {
    for(int jy=1;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++) {
	result(jx, jy, jz) = 2.*(var(jx, jy, jz) - var(jx, jy-1, jz)) / (metric->dy(jx, jy) * sqrt(metric->g_22(jx, jy)) + metric->dy(jx, jy-1) * sqrt(metric->g_22(jx, jy-1)));
      }
    }
  }

  return result;
}

const Field3D Vpar_Grad_par_LCtoC(const Field &v, const Field &f) {
  bindex bx;
  bstencil fval, vval;
  Field3D result;
  
  result.allocate();

  start_index(&bx);
  do {
    f.setStencil(&fval, &bx);
    v.setStencil(&vval, &bx);
    
    // Left side
    result(bx.jx, bx.jy, bx.jz) = (vval.cc >= 0.0) ? vval.cc * fval.ym : vval.cc * fval.cc;
    // Right side
    result(bx.jx, bx.jy, bx.jz) -= (vval.yp >= 0.0) ? vval.yp * fval.cc : vval.yp * fval.yp;
    
  }while(next_index3(&bx));

  return result;
}

const Field3D Grad_par_LtoC(const Field3D &var) {
  Field3D result;
  result.allocate();
  
  Coordinates *metric = mesh->coordinates();
  
  for(int jx=0; jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ngy-1;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++) {
	result(jx, jy, jz) = (var(jx, jy+1, jz) - var(jx, jy, jz)) / (metric->dy(jx, jy) * sqrt(metric->g_22(jx, jy)));
      }
    }
  }
  
  return result;
}

const Field3D Div_par_LtoC(const Field2D &var) {
  Coordinates *metric = mesh->coordinates();
  return metric->Bxy*Grad_par_LtoC(var/metric->Bxy);
}

const Field3D Div_par_LtoC(const Field3D &var) {
  Coordinates *metric = mesh->coordinates();
  return metric->Bxy*Grad_par_LtoC(var/metric->Bxy);
}

const Field3D Div_par_CtoL(const Field2D &var) {
  Coordinates *metric = mesh->coordinates();
  return metric->Bxy*Grad_par_CtoL(var/metric->Bxy);
}

const Field3D Div_par_CtoL(const Field3D &var) {
  Coordinates *metric = mesh->coordinates();
  return metric->Bxy*Grad_par_CtoL(var/metric->Bxy);
}

/*******************************************************************************
 * Grad2_par2
 * second parallel derivative
 *
 * (b dot Grad)(b dot Grad)
 *
 * Note: For parallel Laplacian use LaplacePerp
 *******************************************************************************/

const Field2D Grad2_par2(const Field2D &f) {
  return mesh->coordinates()->Grad2_par2(f);
}

const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc) {
  return mesh->coordinates()->Grad2_par2(f, outloc);
}

/*******************************************************************************
 * Div_par_K_Grad_par
 * Parallel divergence of diffusive flux, K*Grad_par
 *******************************************************************************/

const Field2D Div_par_K_Grad_par(BoutReal kY, Field2D &f) {
  return kY*Grad2_par2(f);
}

const Field3D Div_par_K_Grad_par(BoutReal kY, Field3D &f) {
  return kY*Grad2_par2(f);
}

const Field2D Div_par_K_Grad_par(Field2D &kY, Field2D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(Field2D &kY, Field3D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(Field3D &kY, Field2D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(Field3D &kY, Field3D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

/*******************************************************************************
 * Div_K_perp_Grad_perp
 * Divergence of perpendicular diffusive flux kperp*Grad_perp
 *******************************************************************************/

const Field3D Div_K_perp_Grad_perp(const Field2D &kperp, const Field3D &f) {
  throw BoutException("Div_K_perp_Grad_per not implemented yet");
  Field3D result = 0.0;
  return result;
}

/*******************************************************************************
 * Delp2
 * perpendicular Laplacian operator
 *******************************************************************************/

const Field2D Delp2(const Field2D &f) {
  return mesh->coordinates()->Delp2(f);
}

const Field3D Delp2(const Field3D &f, BoutReal zsmooth) {
  return mesh->coordinates()->Delp2(f);
}

const FieldPerp Delp2(const FieldPerp &f, BoutReal zsmooth) {
  return mesh->coordinates()->Delp2(f);
}

/*******************************************************************************
 * LaplacePerp
 * Full perpendicular Laplacian operator on scalar field
 *
 * Laplace_perp = Laplace - Laplace_par
 *******************************************************************************/

const Field2D Laplace_perp(const Field2D &f) {
  return Laplace(f) - Laplace_par(f);
}

const Field3D Laplace_perp(const Field3D &f) {
  return Laplace(f) - Laplace_par(f);
}

/*******************************************************************************
 * LaplacePar
 * Full parallel Laplacian operator on scalar field
 *
 * LaplacePar(f) = Div( b (b dot Grad(f)) ) 
 *
 *******************************************************************************/

const Field2D Laplace_par(const Field2D &f) {
  return mesh->coordinates()->Laplace_par(f);
}

const Field3D Laplace_par(const Field3D &f) {
  return mesh->coordinates()->Laplace_par(f);
}

/*******************************************************************************
 * Laplacian
 * Full Laplacian operator on scalar field
 *******************************************************************************/

const Field2D Laplace(const Field2D &f) {
  return mesh->coordinates()->Laplace(f);
}

const Field3D Laplace(const Field3D &f) {
  return mesh->coordinates()->Laplace(f);
}

/*******************************************************************************
 * b0xGrad_dot_Grad
 * Terms of form b0 x Grad(phi) dot Grad(A)
 * Used for ExB terms and perturbed B field using A_||
 *******************************************************************************/

const Field2D b0xGrad_dot_Grad(const Field2D &phi, const Field2D &A) {
  Field2D dpdx, dpdy;
  Field2D vx, vy;
  Field2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field2D , Field2D )");
#endif
  
  Coordinates *metric = mesh->coordinates();

  // Calculate phi derivatives
  dpdx = DDX(phi);
  dpdy = DDY(phi);
  
  // Calculate advection velocity
  vx = -metric->g_23*dpdy;
  vy = metric->g_23*dpdx;

  // Upwind A using these velocities
  result = VDDX(vx, A) + VDDY(vy, A);
  result /= metric->J*sqrt(metric->g_22);

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A) {
  Field2D dpdx, dpdy;
  Field2D vx, vy, vz;
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field2D , Field3D )");
#endif

  Coordinates *metric = mesh->coordinates();
  
  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(phi); 
    
    #pragma omp section
    dpdy = DDY(phi);
  }
  
  // Calculate advection velocity
  #pragma omp parallel sections
  {
    #pragma omp section
    vx = -metric->g_23*dpdy;
    
    #pragma omp section
    vy = metric->g_23*dpdx;
    
    #pragma omp section
    vz = metric->g_12*dpdy - metric->g_22*dpdx;
  }

  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += mesh->IntShiftTorsion * vx;
  }

  // Upwind A using these velocities
  
  Field3D ry,rz;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    ry = VDDY(vy, A);

    #pragma omp section
    rz = VDDZ(vz, A);
  }

  result = (result + ry + rz) / (metric->J*sqrt(metric->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &p, const Field2D &A, CELL_LOC outloc) {
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy;
  Field3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field2D )");
#endif

  Coordinates *metric = mesh->coordinates();

  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(p, outloc);
    
    #pragma omp section
    dpdy = DDY(p, outloc);
    
    #pragma omp section
    dpdz = DDZ(p, outloc);
  }

  // Calculate advection velocity
  #pragma omp parallel sections
  {
    #pragma omp section
    vx = metric->g_22*dpdz - metric->g_23*dpdy;
    
    #pragma omp section
    vy = metric->g_23*dpdx - metric->g_12*dpdz;
  }

  // Upwind A using these velocities

  Field3D r2;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    r2 = VDDY(vy, A);
  }

  result = (result + r2) / (metric->J*sqrt(metric->g_22));
  
#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+p.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc) {
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field3D )");
#endif

  Coordinates *metric = mesh->coordinates();

  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(phi, outloc); 
    
    #pragma omp section
    dpdy = DDY(phi, outloc);
    
    #pragma omp section
    dpdz = DDZ(phi, outloc);
  }
  
  // Calculate advection velocity
  #pragma omp parallel sections
  {
    #pragma omp section
    vx = metric->g_22*dpdz - metric->g_23*dpdy;
    
    #pragma omp section
    vy = metric->g_23*dpdx - metric->g_12*dpdz;
    
    #pragma omp section
    vz = metric->g_12*dpdy - metric->g_22*dpdx;
  }

  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += mesh->IntShiftTorsion * vx;
  }

  // Upwind A using these velocities
  
  Field3D ry, rz;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    ry = VDDY(vy, A);
    
    #pragma omp section
    rz = VDDZ(vz, A);
  }
  
  result = (result + ry + rz) / (metric->J*sqrt(metric->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

/*******************************************************************************
 * Poisson bracket
 * Terms of form b0 x Grad(f) dot Grad(g) / B = [f, g]
 *******************************************************************************/

const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method, Solver *solver) {
  Field2D result;
  if( (method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / mesh->coordinates()->Bxy;
  }
  return result;
}

const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method, Solver *solver) {
  Field3D result;
  
  Coordinates *metric = mesh->coordinates();

  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");
    
    // Get current timestep
    BoutReal dt = solver->getCurrentTimestep();

    result.allocate();
    
    int ncz = mesh->ngz - 1;
    for(int x=mesh->xstart;x<=mesh->xend;x++)
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        for(int z=0;z<ncz;z++) {
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;
          
          BoutReal gp, gm;

          // Vx = DDZ(f)
          BoutReal vx = (f(x,y,zp) - f(x,y,zm))/(2.*metric->dz);
          
          // Set stability condition
          solver->setMaxTimestep(metric->dx(x,y) / (fabs(vx) + 1e-16));
          
          // X differencing
          if(vx > 0.0) {
            gp = g(x,y);
            
            gm = g(x-1,y);
            
          }else {
            gp = g(x+1,y);
            
            gm = g(x,y);
          }
          
          result(x,y,z) = vx * (gp - gm) / metric->dx(x,y);
        }
      }
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow. Here as a test

    Field3D fs = f;
    if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
      fs = f.shiftZ(true);
    }
    
    result.allocate();
    int ncz = mesh->ngz - 1;
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<ncz;jz++) {
          int jzp = (jz + 1) % ncz;
          int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = 0.25*( (fs[jx][jy][jzp] - fs[jx][jy][jzm])*
                                (g[jx+1][jy] - g[jx-1][jy]) -
                                (fs[jx+1][jy][jz] - fs[jx-1][jy][jz])*
                                (g[jx][jy] - g[jx][jy]) )
            / (metric->dx[jx][jy] * metric->dz);

          // J+x
          BoutReal Jpx = 0.25*( g[jx+1][jy]*(fs[jx+1][jy][jzp]-fs[jx+1][jy][jzm]) -
                                g[jx-1][jy]*(fs[jx-1][jy][jzp]-fs[jx-1][jy][jzm]) -
                                g[jx][jy]*(fs[jx+1][jy][jzp]-fs[jx-1][jy][jzp]) +
                                g[jx][jy]*(fs[jx+1][jy][jzm]-fs[jx-1][jy][jzm]))
            / (metric->dx[jx][jy] * metric->dz);
          // Jx+
          BoutReal Jxp = 0.25*( g[jx+1][jy]*(fs[jx][jy][jzp]-fs[jx+1][jy][jz]) -
                                g[jx-1][jy]*(fs[jx-1][jy][jz]-fs[jx][jy][jzm]) -
                                g[jx-1][jy]*(fs[jx][jy][jzp]-fs[jx-1][jy][jz]) +
                                g[jx+1][jy]*(fs[jx+1][jy][jz]-fs[jx][jy][jzm]))
            / (metric->dx[jx][jy] * metric->dz);
          
          result[jx][jy][jz] = (Jpp + Jpx + Jxp) / 3.;
        }
    
    if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
      result = result.shiftZ(false); // Shift back
    
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f), g);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / metric->Bxy;
  }
  }
  return result;
}

const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method, Solver *solver) {
  Field3D result;
  switch(method) {
  case BRACKET_CTU:
  case BRACKET_ARAKAWA: 
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(-DDX(f), g);
    break;
  }
  default: {
    // Use full expression with all terms
    Coordinates *metric = mesh->coordinates();
    result = b0xGrad_dot_Grad(f, g) / metric->Bxy;
  }
  }
  return result;
}

const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method, Solver *solver) {

  Coordinates *metric = mesh->coordinates();

  Field3D result;
  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");

    // Get current timestep
    BoutReal dt = solver->getCurrentTimestep();

    result.allocate();
    
    FieldPerp vx, vz;
    vx.allocate();
    vz.allocate();
    
    Field3D fs = f;
    Field3D gs = g;
    if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
      fs = f.shiftZ(true);
      gs = g.shiftZ(true);
    }
    
    int ncz = mesh->ngz - 1;
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      for(int x=1;x<=mesh->ngx-2;x++) {
        for(int z=0;z<ncz;z++) {
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;
          
          // Vx = DDZ(f)
          vx[x][z] = (fs(x,y,zp) - fs(x,y,zm))/(2.*metric->dz);
          // Vz = -DDX(f)
          vz[x][z] = (fs(x-1,y,z) - fs(x+1,y,z))/(0.5*metric->dx(x-1,y) + metric->dx(x,y) + 0.5*metric->dx(x+1,y));
          
          // Set stability condition
          solver->setMaxTimestep(fabs(metric->dx(x,y)) / (fabs(vx(x,z)) + 1e-16));
          solver->setMaxTimestep(metric->dz / (fabs(vz(x,z)) + 1e-16));
        }
      }
      
      // Simplest form: use cell-centered velocities (no divergence included so not flux conservative)
      
      for(int x=mesh->xstart;x<=mesh->xend;x++)
        for(int z=0;z<ncz;z++) {
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;
          
          BoutReal gp, gm;

          // X differencing
          if(vx[x][z] > 0.0) {
            gp = gs(x,y,z)
              + (0.5*dt/metric->dz) * ( (vz(x,z) > 0) ? vz(x,z)*(gs(x,y,zm) - gs(x,y,z)) : vz(x,z)*(gs(x,y,z) - gs(x,y,zp)) );
            
            
            gm = gs(x-1,y,z)
              //+ (0.5*dt/metric->dz) * ( (vz[x-1][z] > 0) ? vz[x-1][z]*(g[x-1][y][zm] - g(x-1,y,z)) : vz[x-1][z]*(g(x-1,y,z) - g[x-1][y][zp]) );
              + (0.5*dt/metric->dz) * ( (vz(x,z) > 0) ? vz(x,z)*(gs(x-1,y,zm) - gs(x-1,y,z)) : vz(x,z)*(gs(x-1,y,z) - gs(x-1,y,zp)) );
            
          }else {
            gp = gs(x+1,y,z)
              //+ (0.5*dt/metric->dz) * ( (vz[x+1][z] > 0) ? vz[x+1][z]*(gs[x+1][y][zm] - gs(x+1,y,z)) : vz[x+1][z]*(gs(x+1,y,z) - gs[x+1][y][zp]) );
              + (0.5*dt/metric->dz) * ( (vz(x,z) > 0) ? vz(x,z)*(gs(x+1,y,zm) - gs(x+1,y,z)) : vz[x][z]*(gs(x+1,y,z) - gs(x+1,y,zp)) );
            
            gm = gs(x,y,z) 
              + (0.5*dt/metric->dz) * ( (vz(x,z) > 0) ? vz(x,z)*(gs(x,y,zm) - gs(x,y,z)) : vz(x,z)*(gs(x,y,z) - gs(x,y,zp)) );
          }
          
          result(x,y,z) = vx(x,z) * (gp - gm) / metric->dx(x,y);
          
          // Z differencing
          if(vz(x,z) > 0.0) {
            gp = gs(x,y,z)
              + (0.5*dt/metric->dx(x,y)) * ( (vx[x][z] > 0) ? vx[x][z]*(gs(x-1,y,z) - gs(x,y,z)) : vx[x][z]*(gs(x,y,z) - gs(x+1,y,z)) );
            
            gm = gs(x,y,zm)
              //+ (0.5*dt/metric->dx(x,y)) * ( (vx[x][zm] > 0) ? vx[x][zm]*(gs[x-1][y][zm] - gs(x,y,zm)) : vx[x][zm]*(gs(x,y,zm) - gs[x+1][y][zm]) );
              + (0.5*dt/metric->dx(x,y)) * ( (vx[x][z] > 0) ? vx[x][z]*(gs[x-1][y][zm] - gs(x,y,zm)) : vx[x][z]*(gs(x,y,zm) - gs(x+1,y,zm)) );
          }else {
            gp = gs(x,y,zp)
              //+ (0.5*dt/metric->dx(x,y)) * ( (vx[x][zp] > 0) ? vx[x][zp]*(gs[x-1][y][zp] - gs[x][y][zp]) : vx[x][zp]*(gs[x][y][zp] - gs[x+1][y][zp]) );
              + (0.5*dt/metric->dx(x,y)) * ( (vx[x][z] > 0) ? vx[x][z]*(gs[x-1][y][zp] - gs[x][y][zp]) : vx[x][z]*(gs(x,y,zp) - gs(x+1,y,zp)) );
            
            gm = gs(x,y,z)
              + (0.5*dt/metric->dx(x,y)) * ( (vx(x,z) > 0) ? vx(x,z)*(gs(x-1,y,z) - gs(x,y,z)) : vx(x,z)*(gs(x,y,z) - gs(x+1,y,z)) );
          }
          
          result(x,y,z) += vz(x,z) * (gp - gm) / metric->dz;
        }
    }
    if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
      result = result.shiftZ(false); // Shift back
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow
    
    result.allocate();
    
    Field3D fs = f;
    Field3D gs = g;
    if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
      fs = f.shiftZ(true);
      gs = g.shiftZ(true);
    }
    
    int ncz = mesh->ngz - 1;
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<ncz;jz++) {
          int jzp = (jz + 1) % ncz;
          int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = 0.25*( (fs[jx][jy][jzp] - fs[jx][jy][jzm])*
                                (gs[jx+1][jy][jz] - gs[jx-1][jy][jz]) -
                                (fs[jx+1][jy][jz] - fs[jx-1][jy][jz])*
                                (gs[jx][jy][jzp] - gs[jx][jy][jzm]) )
            / (metric->dx[jx][jy] * metric->dz);

          // J+x
          BoutReal Jpx = 0.25*( gs[jx+1][jy][jz]*(fs[jx+1][jy][jzp]-fs[jx+1][jy][jzm]) -
                                gs[jx-1][jy][jz]*(fs[jx-1][jy][jzp]-fs[jx-1][jy][jzm]) -
                                gs[jx][jy][jzp]*(fs[jx+1][jy][jzp]-fs[jx-1][jy][jzp]) +
                                gs[jx][jy][jzm]*(fs[jx+1][jy][jzm]-fs[jx-1][jy][jzm]))
            / (metric->dx[jx][jy] * metric->dz);
          // Jx+
          BoutReal Jxp = 0.25*( gs[jx+1][jy][jzp]*(fs[jx][jy][jzp]-fs[jx+1][jy][jz]) -
                                gs[jx-1][jy][jzm]*(fs[jx-1][jy][jz]-fs[jx][jy][jzm]) -
                                gs[jx-1][jy][jzp]*(fs[jx][jy][jzp]-fs[jx-1][jy][jz]) +
                                gs[jx+1][jy][jzm]*(fs[jx+1][jy][jz]-fs[jx][jy][jzm]))
            / (metric->dx(jx,jy) * metric->dz);
          
          result(jx,jy,jz) = (Jpp + Jpx + Jxp) / 3.;
        }
    if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
      result = result.shiftZ(false); // Shift back
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f), g) + VDDZ(-DDX(f), g);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / metric->Bxy;
  }
  }
  return result;
}
