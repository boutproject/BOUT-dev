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
#include <difops.hxx>
#include <vecops.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <fft.hxx>

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

#include <interpolation.hxx>

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


  Field2D result = mesh->Bxy*Grad_par(f/mesh->Bxy);


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


  Field3D result = mesh->Bxy*Grad_par(f/mesh->Bxy, outloc, method);


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
  result.allocate();
  BoutReal ***d = result.getData();

  /*
  bindex bx;
  bstencil f;
  start_index(&bx);
  do {
    var.setStencil(&f, &bx);
    
    d[bx.jx][bx.jy][bx.jz] = (f.cc - f.ym) / dy[bx.jx][bx.jy];
  }while(next_index3(&bx));
  */
  
  // NOTE: Need to calculate one more point than centred vars
  for(int jx=0; jx<mesh->ngx;jx++) {
    for(int jy=1;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++) {
	d[jx][jy][jz] = (var[jx][jy][jz] - var[jx][jy-1][jz]) / (mesh->dy[jx][jy] * sqrt(mesh->g_22[jx][jy]));
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
  
  result.allocate();
  BoutReal ***d = result.getData();

  start_index(&bx);
  do {
    f.setStencil(&fval, &bx);
    v.setStencil(&vval, &bx);
    
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
  
  result.allocate();
  BoutReal ***d = result.getData();

  start_index(&bx);
  do {
    var.setStencil(&f, &bx);
    
    d[bx.jx][bx.jy][bx.jz] = (f.yp - f.cc) / (mesh->dy[bx.jx][bx.jy] * sqrt(mesh->g_22[bx.jx][bx.jy]));
  }while(next_index3(&bx));

  return result;
}

const Field3D Div_par_LtoC(const Field2D &var)
{
  Field3D result = mesh->Bxy*Grad_par_LtoC(var/mesh->Bxy);
  return result;
}

const Field3D Div_par_LtoC(const Field3D &var)
{
  Field3D result = mesh->Bxy*Grad_par_LtoC(var/mesh->Bxy);
  return result;
}

const Field3D Div_par_CtoL(const Field2D &var)
{
  Field3D result = mesh->Bxy*Grad_par_CtoL(var/mesh->Bxy);
  return result;
}

const Field3D Div_par_CtoL(const Field3D &var)
{
  Field3D result = mesh->Bxy*Grad_par_CtoL(var/mesh->Bxy);
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

  Field2D sg;
  Field3D result, r2;
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      sg = sqrt(mesh->g_22);
      sg = DDY(1./sg) / sg;
    }
    
    #pragma omp section
    result = DDY(f);
    
    #pragma omp section
    r2 = D2DY2(f)/mesh->g_22;
  }
  result = sg*result + r2;
  
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

const Field3D Div_par_K_Grad_par(Field3D &kY, Field3D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

/*******************************************************************************
 * Div_K_perp_Grad_perp
 * Divergence of perpendicular diffusive flux kperp*Grad_perp
 *******************************************************************************/

const Field3D Div_K_perp_Grad_perp(const Field2D &kperp, const Field3D &f) {
  
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

const Field3D Delp2(const Field3D &f, BoutReal zsmooth)
{
  Field3D result;
  BoutReal ***fd, ***rd;

#ifdef CHECK
  int msg_pos = msg_stack.push("Delp2( Field3D )");
#endif

  //return mesh->G1*DDX(f) + mesh->G3*DDZ(f) + mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f); //+ 2.0*mesh->g13*D2DXDZ(f)

  // NEW: SOLVE USING FFT

  static dcomplex **ft = (dcomplex**) NULL, **delft;

  result.allocate();

  fd = f.getData();
  rd = result.getData();

  int ncz = mesh->ngz-1;
  
  if(ft == (dcomplex**) NULL) {
    //.allocate memory
    ft = cmatrix(mesh->ngx, ncz/2 + 1);
    delft = cmatrix(mesh->ngx, ncz/2 + 1);
  }
  
  // Loop over all y indices
  for(int jy=0;jy<mesh->ngy;jy++) {

    // Take forward FFT
    
    #pragma omp parallel for
    for(int jx=0;jx<mesh->ngx;jx++)
      ZFFT(fd[jx][jy], mesh->zShift[jx][jy], ft[jx]);

    // Loop over kz
    #pragma omp parallel for
    for(int jz=0;jz<=ncz/2;jz++) {
      BoutReal filter;
      dcomplex a, b, c;
      
      if ((zsmooth > 0.0) && (jz > (int) (zsmooth*((BoutReal) ncz)))) filter=0.0; else filter=1.0;

      // No smoothing in the x direction
      for(int jx=2;jx<(mesh->ngx-2);jx++) {
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
    #pragma omp parallel for
    for(int jx=1;jx<(mesh->ngx-1);jx++) {

      ZFFT_rev(delft[jx], mesh->zShift[jx][jy], rd[jx][jy]);
      rd[jx][jy][ncz] = rd[jx][jy][0];
    }

    // Boundaries
    for(int jz=0;jz<ncz;jz++) {
      rd[0][jy][jz] = 0.0;
      rd[mesh->ngx-1][jy][jz] = 0.0;
    }
  }

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  // Set the output location
  result.setLocation(f.getLocation());

  return result;
}

const FieldPerp Delp2(const FieldPerp &f, BoutReal zsmooth)
{
  FieldPerp result;
  result.allocate();
  
#ifdef CHECK
  int msg_pos = msg_stack.push("Delp2( FieldPerp )");
#endif

  static dcomplex **ft = (dcomplex**) NULL, **delft;
  BoutReal filter;
  
  BoutReal **fd = f.getData();
  BoutReal **rd = result.getData();

  int jy = f.getIndex();
  result.setIndex(jy);
  
  int ncz = mesh->ngz-1;
  
  if(ft == (dcomplex**) NULL) {
    //Allocate memory
    ft = cmatrix(mesh->ngx, ncz/2 + 1);
    delft = cmatrix(mesh->ngx, ncz/2 + 1);
  }
  
  // Take forward FFT
  for(int jx=0;jx<mesh->ngx;jx++)
    ZFFT(fd[jx], mesh->zShift[jx][jy], ft[jx]);

  // Loop over kz
  for(int jz=0;jz<=ncz/2;jz++) {

    if ((zsmooth > 0.0) && (jz > (int) (zsmooth*((BoutReal) ncz)))) filter=0.0; else filter=1.0;
    
    // No smoothing in the x direction
    for(int jx=2;jx<(mesh->ngx-2);jx++) {
      // Perform x derivative
      
      dcomplex a, b, c;
      laplace_tridag_coefs(jx, jy, jz, a, b, c);
      
      delft[jx][jz] = a*ft[jx-1][jz] + b*ft[jx][jz] + c*ft[jx+1][jz];
      delft[jx][jz] *= filter;
    }
  }
  
  // Reverse FFT
  for(int jx=1;jx<(mesh->ngx-1);jx++) {
    ZFFT_rev(delft[jx], mesh->zShift[jx][jy], rd[jx]);
    rd[jx][ncz] = rd[jx][0];
  }

  // Boundaries
  for(int jz=0;jz<ncz;jz++) {
    rd[0][jz] = 0.0;
    rd[mesh->ngx-1][jz] = 0.0;
  }

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

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
  Field2D dpdx, dpdy;
  Field2D vx, vy;
  Field2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field2D , Field2D )");
#endif
  
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
    vx = -mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_23*dpdx;
  }

  // Upwind A using these velocities
  Field2D r2;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    r2 = VDDY(vy, A);
  }
  result += r2;
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
  Field2D dpdx, dpdy;
  Field2D vx, vy, vz;
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field2D , Field3D )");
#endif

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
    vx = -mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_23*dpdx;
    
    #pragma omp section
    vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
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

  result = (result + ry + rz) / (mesh->J*sqrt(mesh->g_22));

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
  Field3D vx, vy;
  Field3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field2D )");
#endif

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
    vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
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

  result = (result + r2) / (mesh->J*sqrt(mesh->g_22));
  
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
    vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
    
    #pragma omp section
    vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
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
  
  result = (result + ry + rz) / (mesh->J*sqrt(mesh->g_22));

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

const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method)
{
  Field2D result;
  if( (method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
  }
  return result;
}

const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method)
{
  Field3D result;
  switch(method) {
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow. Here as a test
    
    result.allocate();
    int ncz = mesh->ngz - 1;
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<ncz;jz++) {
          int jzp = (jz + 1) % ncz;
          int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = 0.25*( (f[jx][jy][jzp] - f[jx][jy][jzm])*
                                (g[jx+1][jy] - g[jx-1][jy]) -
                                (f[jx+1][jy][jz] - f[jx-1][jy][jz])*
                                (g[jx][jy] - g[jx][jy]) )
            / (mesh->dx[jx][jy] * mesh->dz);

          // J+x
          BoutReal Jpx = 0.25*( g[jx+1][jy]*(f[jx+1][jy][jzp]-f[jx+1][jy][jzm]) -
                                g[jx-1][jy]*(f[jx-1][jy][jzp]-f[jx-1][jy][jzm]) -
                                g[jx][jy]*(f[jx+1][jy][jzp]-f[jx-1][jy][jzp]) +
                                g[jx][jy]*(f[jx+1][jy][jzm]-f[jx-1][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          // Jx+
          BoutReal Jxp = 0.25*( g[jx+1][jy]*(f[jx][jy][jzp]-f[jx+1][jy][jz]) -
                                g[jx-1][jy]*(f[jx-1][jy][jz]-f[jx][jy][jzm]) -
                                g[jx-1][jy]*(f[jx][jy][jzp]-f[jx-1][jy][jz]) +
                                g[jx+1][jy]*(f[jx+1][jy][jz]-f[jx][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          
          result[jx][jy][jz] = (Jpp + Jpx + Jxp) / 3.;
        }
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f), g);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
  }
  }
  return result;
}

const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method)
{
  Field3D result;
  switch(method) {
  case BRACKET_ARAKAWA: 
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(-DDX(f), g);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
  }
  }
  return result;
}

const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method)
{
  Field3D result;
  switch(method) {
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow. Here as a test
    
    result.allocate();
    
    int ncz = mesh->ngz - 1;
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++)
        for(int jz=0;jz<ncz;jz++) {
          int jzp = (jz + 1) % ncz;
          int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = 0.25*( (f[jx][jy][jzp] - f[jx][jy][jzm])*
                                (g[jx+1][jy][jz] - g[jx-1][jy][jz]) -
                                (f[jx+1][jy][jz] - f[jx-1][jy][jz])*
                                (g[jx][jy][jzp] - g[jx][jy][jzm]) )
            / (mesh->dx[jx][jy] * mesh->dz);

          // J+x
          BoutReal Jpx = 0.25*( g[jx+1][jy][jz]*(f[jx+1][jy][jzp]-f[jx+1][jy][jzm]) -
                                g[jx-1][jy][jz]*(f[jx-1][jy][jzp]-f[jx-1][jy][jzm]) -
                                g[jx][jy][jzp]*(f[jx+1][jy][jzp]-f[jx-1][jy][jzp]) +
                                g[jx][jy][jzm]*(f[jx+1][jy][jzm]-f[jx-1][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          // Jx+
          BoutReal Jxp = 0.25*( g[jx+1][jy][jzp]*(f[jx][jy][jzp]-f[jx+1][jy][jz]) -
                                g[jx-1][jy][jzm]*(f[jx-1][jy][jz]-f[jx][jy][jzm]) -
                                g[jx-1][jy][jzp]*(f[jx][jy][jzp]-f[jx-1][jy][jz]) +
                                g[jx+1][jy][jzm]*(f[jx+1][jy][jz]-f[jx][jy][jzm]))
            / (mesh->dx[jx][jy] * mesh->dz);
          
          result[jx][jy][jz] = (Jpp + Jpx + Jxp) / 3.;
        }
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f), g) + VDDZ(-DDX(f), g);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
  }
  }
  return result;
}
