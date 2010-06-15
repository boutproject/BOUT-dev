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

#include "globals.h"
#include "difops.h"
#include "utils.h"
#include "derivs.h"
#include "fft.h"

#include "invert_laplace.h" // Delp2 uses same coefficients as inversion code

#include "interpolation.h"

#include <math.h>
#include <stdlib.h>

/*******************************************************************************
 * Grad_par
 * The parallel derivative along unperturbed B-field
 *******************************************************************************/

const Field2D Grad_par(const Field2D &var, CELL_LOC outloc, DIFF_METHOD method)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Grad_par( Field2D )");
#endif


  Field2D result = DDY(var)/sqrt(mesh->g_22); // NOTE: 2D functions not implemented yet


#ifdef TRACK
  result.name = "Grad_par("+var.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field2D Grad_par(const Field2D &var, DIFF_METHOD method, CELL_LOC outloc)
{
  return Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, CELL_LOC outloc, DIFF_METHOD method)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Grad_par( Field3D )");
#endif

  Field3D result;

  result = DDY(var, outloc, method)/sqrt(mesh->g_22);

#ifdef TRACK
  result.name = "Grad_par("+var.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D Grad_par(const Field3D &var, DIFF_METHOD method, CELL_LOC outloc)
{
  return Grad_par(var, outloc, method);
}

/*******************************************************************************
 * Vpar_Grad_par
 * vparallel times the parallel derivative along unperturbed B-field
 *******************************************************************************/

const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f)
{
  return VDDY(v, f)/sqrt(mesh->g_22);
}

const Field3D Vpar_Grad_par(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method)
{
  return VDDY(v, f, outloc, method)/sqrt(mesh->g_22);
}

const Field3D Vpar_Grad_par(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return Vpar_Grad_par(v, f, outloc, method);
}

/*******************************************************************************
 * Div_par
 * parallel divergence operator B \partial_{||} (F/B)
 *******************************************************************************/

const Field2D Div_par(const Field2D &f)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Div_par( Field2D )");
#endif


  Field2D result = Bxy*Grad_par(f/Bxy);


#ifdef TRACK
  result.name = "Div_par("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D Div_par(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Div_par( Field3D )");
#endif


  Field3D result = Bxy*Grad_par(f/Bxy, outloc, method);


#ifdef TRACK
  result.name = "Div_par("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D Div_par(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return Div_par(f, outloc, method);
}

/*******************************************************************************
 * Parallel derivatives converting between left and cell centred
 * NOTE: These are a quick hack to test if this works. The whole staggered grid
 *       thing needs to be thought through.
 *******************************************************************************/

const Field3D Grad_par_CtoL(const Field3D &var)
{
  Field3D result;
  result.Allocate();
  real ***d = result.getData();

  /*
  bindex bx;
  bstencil f;
  start_index(&bx);
  do {
    var.SetStencil(&f, &bx);
    
    d[bx.jx][bx.jy][bx.jz] = (f.cc - f.ym) / dy[bx.jx][bx.jy];
  }while(next_index3(&bx));
  */
  
  // NOTE: Need to calculate one more point than centred vars
  for(int jx=0; jx<ngx;jx++) {
    for(int jy=1;jy<ngy;jy++) {
      for(int jz=0;jz<ngz;jz++) {
	d[jx][jy][jz] = (var[jx][jy][jz] - var[jx][jy-1][jz]) / (dy[jx][jy] * sqrt(mesh->g_22[jx][jy]));
      }
    }
  }

  return result;
}

const Field3D Vpar_Grad_par_LCtoC(const Field &v, const Field &f)
{
  bindex bx;
  bstencil fval, vval;
  Field3D result;
  
  result.Allocate();
  real ***d = result.getData();

  start_index(&bx);
  do {
    f.SetStencil(&fval, &bx);
    v.SetStencil(&vval, &bx);
    
    // Left side
    d[bx.jx][bx.jy][bx.jz] = (vval.cc >= 0.0) ? vval.cc * fval.ym : vval.cc * fval.cc;
    // Right side
    d[bx.jx][bx.jy][bx.jz] -= (vval.yp >= 0.0) ? vval.yp * fval.cc : vval.yp * fval.yp;
    
  }while(next_index3(&bx));

  return result;
}

const Field3D Grad_par_LtoC(const Field &var)
{
  bindex bx;
  bstencil f;
  Field3D result;
  
  result.Allocate();
  real ***d = result.getData();

  start_index(&bx);
  do {
    var.SetStencil(&f, &bx);
    
    d[bx.jx][bx.jy][bx.jz] = (f.yp - f.cc) / (dy[bx.jx][bx.jy] * sqrt(mesh->g_22[bx.jx][bx.jy]));
  }while(next_index3(&bx));

  return result;
}

const Field3D Div_par_LtoC(const Field2D &var)
{
  Field3D result = Bxy*Grad_par_LtoC(var/Bxy);
  return result;
}

const Field3D Div_par_LtoC(const Field3D &var)
{
  Field3D result = Bxy*Grad_par_LtoC(var/Bxy);
  return result;
}

/*******************************************************************************
 * Grad2_par2
 * second parallel derivative
 *******************************************************************************/

const Field2D Grad2_par2(const Field2D &f)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Grad2_par2( Field2D )");
#endif


  Field2D sg = sqrt(mesh->g_22);
  Field2D result = DDY(1./sg)*DDY(f)/sg + D2DY2(f)/mesh->g_22;
  //Field2D result = D2DY2(f)/mesh->g_22;

#ifdef TRACK
  result.name = "Grad2_par2("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D Grad2_par2(const Field3D &f)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Grad2_par2( Field3D )");
#endif

  Field2D sg = sqrt(mesh->g_22);
  Field3D result = DDY(1./sg)*DDY(f)/sg + D2DY2(f)/mesh->g_22;
  
  //Field3D result = D2DY2(f)/mesh->g_22;
#ifdef TRACK
  result.name = "Grad2_par2("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

/*******************************************************************************
 * Div_par_K_Grad_par
 * Parallel divergence of diffusive flux, K*Grad_par
 *******************************************************************************/

const Field2D Div_par_K_Grad_par(Field2D &kY, Field2D &f)
{
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(Field2D &kY, Field3D &f)
{
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(Field3D &kY, Field2D &f)
{
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(Field3D &kY, Field3D &f)
{
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

/*******************************************************************************
 * Delp2
 * perpendicular Laplacian operator
 *******************************************************************************/

const Field2D Delp2(const Field2D &f)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Delp2( Field2D )");
#endif

  Field2D result =  mesh->G1*DDX(f) + mesh->g11*D2DX2(f);

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D Delp2(const Field3D &f, real zsmooth)
{
  Field3D result;
  real ***fd, ***rd;

#ifdef CHECK
  int msg_pos = msg_stack.push("Delp2( Field3D )");
#endif

  //return mesh->G1*DDX(f) + mesh->G3*DDZ(f) + mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f); //+ 2.0*mesh->g13*D2DXDZ(f)

  // NEW: SOLVE USING FFT

  static dcomplex **ft = (dcomplex**) NULL, **delft;
  int jx, jy, jz;
  real filter;
  dcomplex a, b, c;

  result.Allocate();

  fd = f.getData();
  rd = result.getData();

  if(ft == (dcomplex**) NULL) {
    // Allocate memory
    ft = cmatrix(ngx, ncz/2 + 1);
    delft = cmatrix(ngx, ncz/2 + 1);
  }
  
  // Loop over all y indices
  for(jy=0;jy<ngy;jy++) {

    // Take forward FFT
    
    for(jx=0;jx<ngx;jx++)
      ZFFT(fd[jx][jy], zShift[jx][jy], ft[jx]);

    // Loop over kz
    for(jz=0;jz<=ncz/2;jz++) {

      if ((zsmooth > 0.0) && (jz > (int) (zsmooth*((real) ncz)))) filter=0.0; else filter=1.0;

      // No smoothing in the x direction
      for(jx=2;jx<(ngx-2);jx++) {
	// Perform x derivative
	
	laplace_tridag_coefs(jx, jy, jz, a, b, c);

	delft[jx][jz] = a*ft[jx-1][jz] + b*ft[jx][jz] + c*ft[jx+1][jz];
	delft[jx][jz] *= filter;
	
	//Savitzky-Golay 2nd order, 2nd degree in x
        /*
	delft[jx][jz] = coef1*(  0.285714 * (ft[jx-2][jz] + ft[jx+2][jz])
				 - 0.142857 * (ft[jx-1][jz] + ft[jx+1][jz])
				 - 0.285714 * ft[jx][jz] );
	
	delft[jx][jz] -= SQ(kwave)*coef2*ft[jx][jz];
	*/
      }
    }
  
    // Reverse FFT
    for(jx=1;jx<(ngx-1);jx++) {

      ZFFT_rev(delft[jx], zShift[jx][jy], rd[jx][jy]);
      rd[jx][jy][ncz] = rd[jx][jy][0];
    }

    // Boundaries
    for(jz=0;jz<ncz;jz++) {
      rd[0][jy][jz] = 0.0;
      rd[ngx-1][jy][jz] = 0.0;
    }
  }

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  // Set the output location
  result.setLocation(f.getLocation());

  return result;
}

/*******************************************************************************
 * Laplacian
 * Full Laplacian operator
 *******************************************************************************/

const Field2D Laplacian(const Field2D &f)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Laplacian( Field2D )");
#endif

  Field2D result =  mesh->G1*DDX(f) + mesh->G2*DDY(f)
      + mesh->g11*D2DX2(f) + mesh->g22*D2DY2(f);
  // + 2.0*mesh->g12*D2DXDY(f);

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D Laplacian(const Field3D &f)
{
#ifdef CHECK
  int msg_pos = msg_stack.push("Laplacian( Field3D )");
#endif

  Field3D result  = mesh->G1*DDX(f) + mesh->G2*DDY(f) + mesh->G3*DDZ(f)
      + mesh->g11*D2DX2(f) + mesh->g22*D2DY2(f) + mesh->g33*D2DZ2(f);
  // + 2.0*(mesh->g12*D2DXDY(f) + mesh->g13*D2DXDZ(f) + mesh->g23*D2DYDZ(f));

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

/*******************************************************************************
 * b0xGrad_dot_Grad
 * Terms of form b0 x Grad(phi) dot Grad(A)
 * Used for ExB terms and perturbed B field using A_||
 *******************************************************************************/

const Field2D b0xGrad_dot_Grad(const Field2D &phi, const Field2D &A)
{
  Field2D dpdx, dpdy, dpdz;
  Field2D vx, vy, vz;
  Field2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field2D , Field2D )");
#endif
  
  // Calculate phi derivatives
  dpdx = DDX(phi); dpdy = DDY(phi); dpdz = DDZ(phi);
  
  // Calculate advection velocity
  vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
  vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
  vz = mesh->g_12*dpdy - mesh->g_22*dpdx;

  if(ShiftXderivs && IncIntShear) {
    // BOUT-06 style differencing
    vz += IntShiftTorsion * vx;
  }

  // Upwind A using these velocities
  
  result = VDDX(vx, A) + VDDY(vy, A) + VDDZ(vz, A);
  
  result /= mesh->J*sqrt(mesh->g_22);

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A)
{
  Field2D dpdx, dpdy, dpdz;
  Field2D vx, vy, vz;
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field2D , Field3D )");
#endif

  // Calculate phi derivatives
  dpdx = DDX(phi); dpdy = DDY(phi); dpdz = DDZ(phi);
  
  // Calculate advection velocity
  vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
  vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
  vz = mesh->g_12*dpdy - mesh->g_22*dpdx;

  if(ShiftXderivs && IncIntShear) {
    // BOUT-06 style differencing
    vz += IntShiftTorsion * vx;
  }

  // Upwind A using these velocities
  
  result = VDDX(vx, A) + VDDY(vy, A) + VDDZ(vz, A);

  result /= mesh->J*sqrt(mesh->g_22);

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &p, const Field2D &A, CELL_LOC outloc)
{
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field2D )");
#endif

  // Calculate phi derivatives
  dpdx = DDX(p, outloc);
  dpdy = DDY(p, outloc);
  dpdz = DDZ(p, outloc);

  // Calculate advection velocity
  vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
  vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
  vz = mesh->g_12*dpdy - mesh->g_22*dpdx;

  if(ShiftXderivs && IncIntShear) {
    // BOUT-06 style differencing
    vz += IntShiftTorsion * vx;
  }

  // Upwind A using these velocities

  result = VDDX(vx, A) + VDDY(vy, A) + VDDZ(vz, A);

  result /= mesh->J*sqrt(mesh->g_22);
  
#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+p.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc)
{
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field3D )");
#endif

  // Calculate phi derivatives
  dpdx = DDX(phi, outloc); dpdy = DDY(phi, outloc); dpdz = DDZ(phi, outloc);
  
  // Calculate advection velocity
  vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
  vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
  vz = mesh->g_12*dpdy - mesh->g_22*dpdx;

  if(ShiftXderivs && IncIntShear) {
    // BOUT-06 style differencing
    vz += IntShiftTorsion * vx;
  }

  // Upwind A using these velocities
  
  result = VDDX(vx, A) + VDDY(vy, A) + VDDZ(vz, A);
  
  result /= mesh->J*sqrt(mesh->g_22);

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}
