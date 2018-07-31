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
#include <unused.hxx>

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
  ASSERT1(apar.getMesh() == f.getMesh());

  Mesh *mesh = apar.getMesh();

  Field3D result(mesh);
  result.allocate();
  
  int ncz = mesh->LocalNz;

  Coordinates *metric = mesh->coordinates();

  Field3D gys(mesh);
  gys.allocate();

  // Need Y derivative everywhere
  for(int x=1;x<=mesh->LocalNx-2;x++)
    for(int y=1;y<=mesh->LocalNy-2;y++)
      for(int z=0;z<ncz;z++) {
        gys(x, y, z) = (f.yup()(x, y+1, z) - f.ydown()(x, y-1, z))/(0.5*metric->dy(x, y+1) + metric->dy(x, y) + 0.5*metric->dy(x, y-1));
      }
  
  for(int x=1;x<=mesh->LocalNx-2;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      BoutReal by = 1./sqrt(metric->g_22(x, y));
      for(int z=0;z<ncz;z++) {
        // Z indices zm and zp
        int zm = (z - 1 + ncz) % ncz;
        int zp = (z + 1) % ncz;
        
        // bx = -DDZ(apar)
        BoutReal bx = (apar(x, y, zm) - apar(x, y, zp))/(2.*metric->dz);
        // bz = DDX(f)
        BoutReal bz = (apar(x+1, y, z) - apar(x-1, y, z))/(0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y));
        
	// Now calculate (bx*d/dx + by*d/dy + bz*d/dz) f
        
        // Length dl for predictor
        BoutReal dl = fabs(metric->dx(x, y)) / (fabs(bx) + 1e-16);
        dl = BOUTMIN(dl, fabs(metric->dy(x, y)) / (fabs(by) + 1e-16));
        dl = BOUTMIN(dl, metric->dz / (fabs(bz) + 1e-16));
        
	BoutReal fp, fm;
        
        // X differencing
        fp = f(x+1, y, z)
          + (0.25*dl/metric->dz) * bz * (f(x+1, y, zm) - f(x+1, y, zp))
          - 0.5*dl * by * gys(x+1, y, z);
        
        fm = f(x-1, y, z)
          + (0.25*dl/metric->dz) * bz * (f(x-1, y, zm) - f(x-1, y, zp))
          - 0.5*dl * by * gys(x-1, y, z);
        
        result(x, y, z) = bx * (fp - fm) / (0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y));

	// Z differencing
        
        fp = f(x, y, zp)
          + (0.25*dl/metric->dx(x, y)) * bx * (f(x-1, y, zp) - f(x+1, y, zp))
          - 0.5*dl * by * gys(x, y, zp);
        
        fm = f(x, y, zm)
          + (0.25*dl/metric->dx(x, y)) * bx * (f(x-1,y,zm) - f(x+1, y, zm))
          - 0.5*dl * by * gys(x, y, zm);

        result(x, y, z) += bz * (fp - fm) / (2.*metric->dz);

        // Y differencing
        
        fp = f.yup()(x,y+1,z)
          - 0.5*dl * bx * (f.yup()(x+1, y+1, z) - f.yup()(x-1, y+1, z))/(0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y))
          
          + (0.25*dl/metric->dz) * bz * (f.yup()(x,y+1,zm) - f.yup()(x,y+1,zp));
        
        fm = f.ydown()(x,y-1,z)
          - 0.5*dl * bx * (f.ydown()(x+1, y-1, z) - f.ydown()(x-1, y-1, z))/(0.5*metric->dx(x-1, y) + metric->dx(x, y) + 0.5*metric->dx(x+1, y))
          + (0.25*dl/metric->dz) * bz * (f.ydown()(x,y-1,zm) - f.ydown()(x,y-1,zp));

        result(x,y,z) += by * (fp - fm) / (0.5*metric->dy(x,y-1) + metric->dy(x,y) + 0.5*metric->dy(x,y+1));
      }
    }
  }
  
  ASSERT2(result.getLocation() == f.getLocation());

  return result;
}

/*******************************************************************************
* Vpar_Grad_par
* vparallel times the parallel derivative along unperturbed B-field
*******************************************************************************/

const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f) {
  return mesh->coordinates()->Vpar_Grad_par(v, f);
}

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return mesh->coordinates()->Vpar_Grad_par(v, f, outloc, method);
}

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
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

const Field3D Div_par(const Field3D &f, const Field3D &v) {
  ASSERT1(f.getMesh() == v.getMesh());
  ASSERT2(v.getLocation() == f.getLocation());

  // Parallel divergence, using velocities at cell boundaries
  // Note: Not guaranteed to be flux conservative
  Mesh *mesh = f.getMesh();
  CELL_LOC loc = f.getLocation();

  Field3D result(mesh);
  result.allocate();

  Coordinates *coord = mesh->coordinates();
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      for(int k=0;k<mesh->LocalNz;k++) {
	
	// Value of f and v at left cell face
	BoutReal fL = 0.5*(f(i,j,k) + f.ydown()(i,j-1,k));
	BoutReal vL = 0.5*(v(i,j,k) + v.ydown()(i,j-1,k));
	
	BoutReal fR = 0.5*(f(i,j,k) + f.yup()(i,j+1,k));
	BoutReal vR = 0.5*(v(i,j,k) + v.yup()(i,j+1,k));
	
        // Calculate flux at right boundary (y+1/2)
	BoutReal fluxRight = fR * vR * (coord->J.get(loc).operator()(i,j) + coord->J.get(loc).operator()(i,j+1)) / (sqrt(coord->g_22.get(loc).operator()(i,j))+ sqrt(coord->g_22.get(loc).operator()(i,j+1)));
	
        // Calculate at left boundary (y-1/2)
	BoutReal fluxLeft = fL * vL * (coord->J.get(loc).operator()(i,j) + coord->J.get(loc).operator()(i,j-1)) / (sqrt(coord->g_22.get(loc).operator()(i,j)) + sqrt(coord->g_22.get(loc).operator()(i,j-1)));
	
	result(i,j,k)   = (fluxRight - fluxLeft) / (coord->dy.get(loc).operator()(i,j)*coord->J(i,j));
      }
    }

  result.setLocation(f.getLocation());

  return result;
}

//////// Flux methods

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  Coordinates *metric = mesh->coordinates();
  Field3D result = -metric->Bxy*FDDY(v, f/metric->Bxy, outloc, method);
  result /= sqrt(metric->g_22.get(result.getLocation()));
  return result;
}

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
  return Div_par_flux(v,f, outloc, method);
}

/*******************************************************************************
* Parallel derivatives converting between left and cell centred
* NOTE: These are a quick hack to test if this works. The whole staggered grid
*       thing needs to be thought through.
*******************************************************************************/

const Field3D Grad_par_CtoL(const Field3D &var) {
  ASSERT1(var.getLocation() == CELL_CENTRE);

  Mesh *mesh = var.getMesh();
  Field3D result(mesh);
  result.allocate();

  Coordinates *metric = mesh->coordinates();

  if (var.hasYupYdown()) {
    // NOTE: Need to calculate one more point than centred vars
    for(int jx=0; jx<mesh->LocalNx;jx++) {
      for(int jy=1;jy<mesh->LocalNy;jy++) {
        for(int jz=0;jz<mesh->LocalNz;jz++) {
          result(jx, jy, jz) = 2.*(var(jx, jy, jz) - var.ydown()(jx, jy-1, jz)) / (metric->dy.get(CELL_YLOW).operator()(jx, jy) * sqrt(metric->g_22.get(CELL_YLOW).operator()(jx, jy)) + metric->dy.get(CELL_YLOW).operator()(jx, jy-1) * sqrt(metric->g_22.get(CELL_YLOW).operator()(jx, jy-1)));
      }
      }
    }
  } else {
    // No yup/ydown fields, so transform to cell centred
    Field3D var_fa = mesh->toFieldAligned(var);
    
    for(int jx=0; jx<mesh->LocalNx;jx++) {
      for(int jy=1;jy<mesh->LocalNy;jy++) {
        for(int jz=0;jz<mesh->LocalNz;jz++) {
          result(jx, jy, jz) = 2.*(var_fa(jx, jy, jz) - var_fa(jx, jy-1, jz)) / (metric->dy.get(CELL_YLOW).operator()(jx, jy) * sqrt(metric->g_22.get(CELL_YLOW).operator()(jx, jy)) + metric->dy.get(CELL_YLOW).operator()(jx, jy-1) * sqrt(metric->g_22.get(CELL_YLOW).operator()(jx, jy-1)));
        }
      }
    }

    result = mesh->fromFieldAligned(result);
  }

  result.setLocation(CELL_YLOW);
  return result;
}

const Field2D Grad_par_CtoL(const Field2D &var) {
  Field2D result;
  result.allocate();
  
  Coordinates *metric = mesh->coordinates();

  for(auto &i : result.region(RGN_NOBNDRY)) {
    result[i] = (var[i] - var[i.ym()]) / (metric->dy[i] * sqrt(metric->g_22[i]));
  }
  
  return result;
}


const Field3D Vpar_Grad_par_LCtoC(const Field3D &v, const Field3D &f, REGION region) {
  ASSERT1(v.getMesh() == f.getMesh());
  ASSERT1(v.getLocation() == CELL_YLOW);
  ASSERT1(f.getLocation() == CELL_CENTRE);

  stencil fval, vval;
  Field3D result(v.getMesh());

  result.allocate();

  bool vUseUpDown = (v.hasYupYdown() && ((&v.yup() != &v) || (&v.ydown() != &v)));
  bool fUseUpDown = (f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f)));

  if (vUseUpDown && fUseUpDown) {
    // Both v and f have up/down fields

    fval.mm = nan("");
    fval.pp = nan("");
    vval.mm = nan("");
    vval.pp = nan("");
    for (const auto &i : result.region(region)) {

      vval.m = v.ydown()[i.ym()];
      vval.c = v[i];
      vval.p = v.yup()[i.yp()];

      fval.m = f.ydown()[i.ym()];
      fval.c = f[i];
      fval.p = f.yup()[i.yp()];

      // Left side
      result[i] = (vval.c >= 0.0) ? vval.c * fval.m : vval.c * fval.c;
      // Right side
      result[i] -= (vval.p >= 0.0) ? vval.p * fval.c : vval.p * fval.p;
    }
  }
  else if (vUseUpDown) {
    // Only v has up/down fields
    // f must shift to field aligned coordinates
    Field3D f_fa = mesh->toFieldAligned(f);

    vval.mm = nan("");
    vval.pp = nan("");

    for (const auto &i : result.region(region)) {

      fval.mm = f_fa[i.offset(0, -2, 0)];
      fval.m = f_fa[i.ym()];
      fval.c = f_fa[i];
      fval.p = f_fa[i.yp()];
      fval.pp = f_fa[i.offset(0, 2, 0)];

      vval.m = v.ydown()[i.ym()];
      vval.c = v[i];
      vval.p = v.yup()[i.yp()];

      // Left side
      result[i] = (vval.c >= 0.0) ? vval.c * fval.m : vval.c * fval.c;
      // Right side
      result[i] -= (vval.p >= 0.0) ? vval.p * fval.c : vval.p * fval.p;
    }
  }
  else if (fUseUpDown) {
    // Only f has up/down fields
    // v must shift to field aligned coordinates
    Field3D v_fa = mesh->toFieldAligned(v);

    stencil vval;

    stencil fval;
    fval.mm = nan("");
    fval.pp = nan("");

    for (const auto &i : result.region(region)) {

      fval.m = f.ydown()[i.ym()];
      fval.c = f[i];
      fval.p = f.yup()[i.yp()];

      vval.mm = v_fa[i.offset(0,-2,0)];
      vval.m = v_fa[i.ym()];
      vval.c = v_fa[i];
      vval.p = v_fa[i.yp()];
      vval.pp = v_fa[i.offset(0,2,0)];

      // Left side
      result[i] = (vval.c >= 0.0) ? vval.c * fval.m : vval.c * fval.c;
      // Right side
      result[i] -= (vval.p >= 0.0) ? vval.p * fval.c : vval.p * fval.p;
    }
  }
  else {
    // Both must shift to field aligned
    Field3D v_fa = mesh->toFieldAligned(v);
    Field3D f_fa = mesh->toFieldAligned(f);

    for (const auto &i : result.region(region)) {

      fval.mm = f_fa[i.offset(0,-2,0)];
      fval.m = f_fa[i.ym()];
      fval.c = f_fa[i];
      fval.p = f_fa[i.yp()];
      fval.pp = f_fa[i.offset(0,2,0)];

      vval.mm = v_fa[i.offset(0,-2,0)];
      vval.m = v_fa[i.ym()];
      vval.c = v_fa[i];
      vval.p = v_fa[i.yp()];
      vval.pp = v_fa[i.offset(0,2,0)];

      // Left side
      result[i] = (vval.c >= 0.0) ? vval.c * fval.m : vval.c * fval.c;
      // Right side
      result[i] -= (vval.p >= 0.0) ? vval.p * fval.c : vval.p * fval.p;
    }
  }

  result.setLocation(CELL_CENTRE);
  return result;
}

const Field3D Grad_par_LtoC(const Field3D &var) {
  if (mesh->StaggerGrids) {
    ASSERT1(var.getLocation() == CELL_YLOW);
  } 

  Field3D result(var.getMesh());
  result.allocate();

  Coordinates *metric = var.getMesh()->coordinates();

  if (var.hasYupYdown()) {
    for (auto &i : result.region(RGN_NOBNDRY)) {
      result[i] = (var.yup()[i.yp()] - var[i]) / (metric->dy[i]*sqrt(metric->g_22.get(CELL_CENTRE).operator[](i)));
    }
  } else {
    // No yup/ydown field, so transform to field aligned

    Field3D var_fa = var.getMesh()->toFieldAligned(var);

    for(auto &i : result.region(RGN_NOBNDRY)) {
      result[i] = (var_fa[i.yp()] - var_fa[i]) / (metric->dy[i]*sqrt(metric->g_22.get(CELL_CENTRE).operator[](i)));
    }
    result = var.getMesh()->fromFieldAligned(result);
  }

  result.setLocation(CELL_CENTRE);
  return result;
}

const Field2D Grad_par_LtoC(const Field2D &var) {
  Field2D result;
  result.allocate();
  
  Coordinates *metric = mesh->coordinates();

  for(auto &i : result.region(RGN_NOBNDRY)) {
    result[i] = (var[i.yp()] - var[i]) / (metric->dy[i] * sqrt(metric->g_22[i]));
  }
  
  return result;
}

const Field2D Div_par_LtoC(const Field2D &var) {
  Coordinates *metric = mesh->coordinates();
  return metric->Bxy*Grad_par_LtoC(var/metric->Bxy);
}

const Field3D Div_par_LtoC(const Field3D &var) {
  Field3D result;
  result.allocate();

  Coordinates *metric = mesh->coordinates();

  // NOTE: Need to calculate one more point than centred vars
  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 0; jy < mesh->LocalNy - 1; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        result(jx, jy, jz) = metric->Bxy(jx, jy) * 2. *
                             (var.yup()(jx, jy + 1, jz) / metric->Bxy(jx, jy + 1) -
                              var(jx, jy, jz) / metric->Bxy(jx, jy)) /
                             (metric->dy(jx, jy) * sqrt(metric->g_22(jx, jy)) +
                              metric->dy(jx, jy - 1) * sqrt(metric->g_22(jx, jy - 1)));
      }
    }
  }

  return result;
}

const Field2D Div_par_CtoL(const Field2D &var) {
  Coordinates *metric = mesh->coordinates();
  return metric->Bxy * Grad_par_CtoL(var / metric->Bxy);
}

const Field3D Div_par_CtoL(const Field3D &var) {
  Field3D result;
  result.allocate();

  Coordinates *metric = mesh->coordinates();

  // NOTE: Need to calculate one more point than centred vars
  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 1; jy < mesh->LocalNy; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        result(jx, jy, jz) = metric->Bxy(jx, jy) * 2. *
                             (var(jx, jy, jz) / metric->Bxy(jx, jy) -
                              var.ydown()(jx, jy - 1, jz) / metric->Bxy(jx, jy - 1)) /
                             (metric->dy(jx, jy) * sqrt(metric->g_22(jx, jy)) +
                              metric->dy(jx, jy - 1) * sqrt(metric->g_22(jx, jy - 1)));
      }
    }
  }

  return result;
}

/*******************************************************************************
* Grad2_par2
* second parallel derivative
*
* (b dot Grad)(b dot Grad)
*
* Note: For parallel Laplacian use LaplacePar
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

const Field2D Div_par_K_Grad_par(BoutReal kY, const Field2D &f) {
  return kY*Grad2_par2(f);
}

const Field3D Div_par_K_Grad_par(BoutReal kY, const Field3D &f) {
  return kY*Grad2_par2(f);
}

const Field2D Div_par_K_Grad_par(const Field2D &kY, const Field2D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(const Field2D &kY, const Field3D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field2D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field3D &f) {
  return kY*Grad2_par2(f) + Div_par(kY)*Grad_par(f);
}

/*******************************************************************************
* Delp2
* perpendicular Laplacian operator
*******************************************************************************/

const Field2D Delp2(const Field2D &f) {
  return mesh->coordinates()->Delp2(f);
}

const Field3D Delp2(const Field3D &f, BoutReal UNUSED(zsmooth)) {
  return mesh->coordinates()->Delp2(f);
}

const FieldPerp Delp2(const FieldPerp &f, BoutReal UNUSED(zsmooth)) {
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
  ASSERT1(phi.getMesh() == A.getMesh());
  ASSERT1(phi.getLocation() == A.getLocation());
  
  TRACE("b0xGrad_dot_Grad( Field2D , Field2D )");
  
  ASSERT1(phi.getMesh() == A.getMesh());

  Mesh * mesh = phi.getMesh();
  Coordinates *metric = mesh->coordinates();
  CELL_LOC loc = phi.getLocation();

  // Calculate phi derivatives
  Field2D dpdx = DDX(phi);
  Field2D dpdy = DDY(phi);
  
  // Calculate advection velocity
  Field2D vx = -metric->g_23*dpdy;
  Field2D vy = metric->g_23*dpdx;

  // Upwind A using these velocities
  Field2D result = VDDX(vx, A) + VDDY(vy, A);
  result /= metric->J*sqrt(metric->g_22.get(loc));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A) {
  TRACE("b0xGrad_dot_Grad( Field2D , Field3D )");

  ASSERT1(phi.getMesh() == A.getMesh());
  ASSERT1(phi.getLocation() == A.getLocation())

  Mesh *mesh = phi.getMesh();

  Coordinates *metric = mesh->coordinates();
  CELL_LOC loc = A.getLocation();
  
  // Calculate phi derivatives
  Field2D dpdx = DDX(phi);
  Field2D dpdy = DDY(phi);

  // Calculate advection velocity
  Field2D vx = -metric->g_23 * dpdy;
  Field2D vy = metric->g_23 * dpdx;
  Field2D vz = metric->g_12 * dpdy - metric->g_22 * dpdx;

  if(mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion * vx;
  }

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A) + VDDY(vy, A) + VDDZ(vz, A);

  result /= (metric->J*sqrt(metric->g_22.get(loc)));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif

  ASSERT2(result.getLocation() == A.getLocation());
  
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &p, const Field2D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field3D , Field2D )");

  ASSERT1(p.getMesh() == A.getMesh());

  Coordinates *metric = p.getMesh()->coordinates();

  // Calculate phi derivatives
  Field3D dpdx = DDX(p, outloc);
  Field3D dpdy = DDY(p, outloc);
  Field3D dpdz = DDZ(p, outloc);

  // Calculate advection velocity
  Field3D vx = metric->g_22 * dpdz - metric->g_23 * dpdy;
  Field3D vy = metric->g_23 * dpdx - metric->g_12 * dpdz;

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);

  result /=  (metric->J*sqrt(metric->g_22.get(result.getLocation())));
  
#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+p.name+","+A.name+")";
#endif
  
  ASSERT2(result.getLocation() == A.getLocation());

  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field3D , Field3D )");

  ASSERT1(phi.getMesh() == A.getMesh());

  Mesh * mesh = phi.getMesh();

  Coordinates *metric = mesh->coordinates();

  // Calculate phi derivatives
  Field3D dpdx = DDX(phi, outloc);
  Field3D dpdy = DDY(phi, outloc);
  Field3D dpdz = DDZ(phi, outloc);

  // Calculate advection velocity
  Field3D vx = metric->g_22 * dpdz - metric->g_23 * dpdy;
  Field3D vy = metric->g_23 * dpdx - metric->g_12 * dpdz;
  Field3D vz = metric->g_12 * dpdy - metric->g_22 * dpdx;

  if(mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion * vx;
  }

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc) + VDDZ(vz, A, outloc);

  result /=  (metric->J*sqrt(metric->g_22.get(result.getLocation())));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == A.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

/*******************************************************************************
 * Poisson bracket
 * Terms of form b0 x Grad(f) dot Grad(g) / B = [f, g]
 *******************************************************************************/

/*!
 * Calculate location of result
 */
CELL_LOC bracket_location(const CELL_LOC &f_loc, const CELL_LOC &g_loc, const CELL_LOC &outloc) {
  if(!mesh->StaggerGrids)
    return CELL_CENTRE;

  if(outloc == CELL_DEFAULT){
    // Check that f and g are in the same location
    if (f_loc != g_loc){
      throw BoutException("Bracket currently requires both fields to have the same cell location");
    }else {
      return f_loc;      	  // Location of result
    }
  }
  
  // Check that f, and g are in the same location as the specified output location
  if(f_loc != g_loc || f_loc != outloc){
    throw BoutException("Bracket currently requires the location of both fields and the output locaton to be the same");
  }
  
  return outloc;      	  // Location of result
}

const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *UNUSED(solver)) {
  TRACE("bracket(Field2D, Field2D)");

  ASSERT1(f.getMesh() == g.getMesh());

  Field2D result(f.getMesh());

  // Sort out cell locations
  CELL_LOC result_loc = bracket_location(f.getLocation(), g.getLocation(), outloc);
  
  if( (method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g) / mesh->coordinates()->Bxy;
  }
  result.setLocation(result_loc);
  return result;
}

const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("bracket(Field3D, Field2D)");

  ASSERT1(f.getMesh() == g.getMesh());

  Mesh *mesh = f.getMesh();

  Field3D result(mesh);

  Coordinates *metric = mesh->coordinates();

  CELL_LOC result_loc = bracket_location(f.getLocation(), g.getLocation(), outloc);
  
  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");
    
    result.allocate();
    
    int ncz = mesh->LocalNz;
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

    result.allocate();
    const int ncz = mesh->LocalNz;
    const BoutReal partialFactor = 1.0/(12 * metric->dz);
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++){
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++){
	const BoutReal spacingFactor = partialFactor / metric->dx(jx,jy);
	for(int jz=0;jz<ncz;jz++) {
	  const int jzp = jz+1 < ncz ? jz + 1 : 0;
	  //Above is alternative to const int jzp = (jz + 1) % ncz;
	  const int jzm = jz-1 >=  0 ? jz - 1 : ncz-1;
	  //Above is alternative to const int jzmTmp = (jz - 1 + ncz) % ncz;

          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = ( (f(jx,jy,jzp) - f(jx,jy,jzm))*
			   (g(jx+1,jy) - g(jx-1,jy)) -
			   (f(jx+1,jy,jz) - f(jx-1,jy,jz))*
			   (g(jx,jy) - g(jx,jy)) );
      
          // J+x
          BoutReal Jpx = ( g(jx+1,jy)*(f(jx+1,jy,jzp)-f(jx+1,jy,jzm)) -
			   g(jx-1,jy)*(f(jx-1,jy,jzp)-f(jx-1,jy,jzm)) -
			   g(jx,jy)*(f(jx+1,jy,jzp)-f(jx-1,jy,jzp)) +
			   g(jx,jy)*(f(jx+1,jy,jzm)-f(jx-1,jy,jzm)));

          // Jx+
          BoutReal Jxp = ( g(jx+1,jy)*(f(jx,jy,jzp)-f(jx+1,jy,jz)) -
			   g(jx-1,jy)*(f(jx-1,jy,jz)-f(jx,jy,jzm)) -
			   g(jx-1,jy)*(f(jx,jy,jzp)-f(jx-1,jy,jz)) +
			   g(jx+1,jy)*(f(jx+1,jy,jz)-f(jx,jy,jzm)));
          
          result(jx,jy,jz) = (Jpp + Jpx + Jxp) * spacingFactor;
	  }
	}
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
    result = b0xGrad_dot_Grad(f, g) / metric->Bxy;
  }
  }
  result.setLocation(result_loc);
  return result;
}

const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("bracket(Field2D, Field3D)");

  ASSERT1(f.getMesh() == g.getMesh());

  Mesh *mesh = f.getMesh();

  Field3D result(mesh);

  CELL_LOC result_loc = bracket_location(f.getLocation(), g.getLocation(), outloc);

  switch(method) {
  case BRACKET_CTU:
    throw BoutException("Bracket method CTU is not yet implemented for [2d,3d] fields.");
    break;
  case BRACKET_ARAKAWA: 
    // It is symmetric, therefore we can return -[3d,2d]
    return -bracket(g,f,method,outloc,solver);
    break;
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
  result.setLocation(result_loc) ;
  
  return result;
}

const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("Field3D, Field3D");

  ASSERT1(f.getMesh() == g.getMesh());

  Mesh *mesh = f.getMesh();

  Coordinates *metric = mesh->coordinates();

  Field3D result(mesh);

  CELL_LOC result_loc = bracket_location(f.getLocation(), g.getLocation(), outloc);

  if (mesh->GlobalNx == 1 || mesh->GlobalNz == 1) {
    result=0;
    result.setLocation(result_loc);
    return result;
  }
  
  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");

    // Get current timestep
    BoutReal dt = solver->getCurrentTimestep();
    
    result.allocate();

    FieldPerp vx(mesh), vz(mesh);
    vx.allocate();
    vz.allocate();
    
    int ncz = mesh->LocalNz;
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      for(int x=1;x<=mesh->LocalNx-2;x++) {
        for(int z=0;z<ncz;z++) {
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;
          
          // Vx = DDZ(f)
          vx(x,z) = (f(x,y,zp) - f(x,y,zm))/(2.*metric->dz);
          // Vz = -DDX(f)
          vz(x,z) = (f(x-1,y,z) - f(x+1,y,z))/(0.5*metric->dx(x-1,y) + metric->dx(x,y) + 0.5*metric->dx(x+1,y));
          
          // Set stability condition
          solver->setMaxTimestep(fabs(metric->dx(x,y)) / (fabs(vx(x,z)) + 1e-16));
          solver->setMaxTimestep(metric->dz / (fabs(vz(x,z)) + 1e-16));
        }
      }
      
      // Simplest form: use cell-centered velocities (no divergence included so not flux conservative)
      
      for(int x=mesh->xstart;x<=mesh->xend;x++)
        for (int z = 0; z < ncz; z++) {
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;

          BoutReal gp, gm;

          // X differencing
          if (vx(x, z) > 0.0) {
            gp = g(x, y, z) +
                 (0.5 * dt / metric->dz) * ((vz(x, z) > 0)
                                                ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                                : vz(x, z) * (g(x, y, z) - g(x, y, zp)));

            gm = g(x - 1, y, z) +
                 (0.5 * dt / metric->dz) *
                     ((vz(x, z) > 0) ? vz(x, z) * (g(x - 1, y, zm) - g(x - 1, y, z))
                                     : vz(x, z) * (g(x - 1, y, z) - g(x - 1, y, zp)));

          } else {
            gp = g(x + 1, y, z) +
                 (0.5 * dt / metric->dz) *
                     ((vz(x, z) > 0) ? vz(x, z) * (g(x + 1, y, zm) - g(x + 1, y, z))
                                     : vz[x][z] * (g(x + 1, y, z) - g(x + 1, y, zp)));

            gm = g(x, y, z) +
                 (0.5 * dt / metric->dz) * ((vz(x, z) > 0)
                                                ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                                : vz(x, z) * (g(x, y, z) - g(x, y, zp)));
          }

          result(x, y, z) = vx(x, z) * (gp - gm) / metric->dx(x, y);

          // Z differencing
          if (vz(x, z) > 0.0) {
            gp = g(x, y, z) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx[x][z] > 0) ? vx[x][z] * (g(x - 1, y, z) - g(x, y, z))
                                     : vx[x][z] * (g(x, y, z) - g(x + 1, y, z)));

            gm = g(x, y, zm) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, zm) - g(x, y, zm))
                                     : vx(x, z) * (g(x, y, zm) - g(x + 1, y, zm)));
          } else {
            gp = g(x, y, zp) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, zp) - g(x, y, zp))
                                     : vx(x, z) * (g(x, y, zp) - g(x + 1, y, zp)));

            gm = g(x, y, z) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, z) - g(x, y, z))
                                     : vx(x, z) * (g(x, y, z) - g(x + 1, y, z)));
          }

          result(x, y, z) += vz(x, z) * (gp - gm) / metric->dz;
        }
    }
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow
    
    result.allocate();
    
    const int ncz = mesh->LocalNz;
    const BoutReal partialFactor = 1.0/(12 * metric->dz);

    // We need to discard const qualifier in order to manipulate
    // storage array directly
    Field3D f_temp = f;
    Field3D g_temp = g;

    for(int jx=mesh->xstart;jx<=mesh->xend;jx++){
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++){
        const BoutReal spacingFactor = partialFactor / metric->dx(jx, jy);
        const BoutReal *Fxm = f_temp(jx-1, jy);
        const BoutReal *Fx  = f_temp(jx,   jy);
        const BoutReal *Fxp = f_temp(jx+1, jy);
        const BoutReal *Gxm = g_temp(jx-1, jy);
        const BoutReal *Gx  = g_temp(jx,   jy);
        const BoutReal *Gxp = g_temp(jx+1, jy);
        for(int jz=0;jz<ncz;jz++) {
	  const int jzp = jz+1 < ncz ? jz + 1 : 0;
	  //Above is alternative to const int jzp = (jz + 1) % ncz;
	  const int jzm = jz-1 >=  0 ? jz - 1 : ncz-1;
	  //Above is alternative to const int jzm = (jz - 1 + ncz) % ncz;

          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = ((Fx[jzp] - Fx[jzm])*(Gxp[jz] - Gxm[jz]) - 
			  (Fxp[jz] - Fxm[jz])*(Gx[jzp] - Gx[jzm]));

          // J+x
          BoutReal Jpx = ( Gxp[jz]*(Fxp[jzp]-Fxp[jzm]) -
			   Gxm[jz]*(Fxm[jzp]-Fxm[jzm]) -
			   Gx[jzp]*(Fxp[jzp]-Fxm[jzp]) +
			   Gx[jzm]*(Fxp[jzm]-Fxm[jzm])) ;

          // Jx+
          BoutReal Jxp = ( Gxp[jzp]*(Fx[jzp]-Fxp[jz]) -
			   Gxm[jzm]*(Fxm[jz]-Fx[jzm]) -
			   Gxm[jzp]*(Fx[jzp]-Fxm[jz]) +
			   Gxp[jzm]*(Fxp[jz]-Fx[jzm]));
			  
          result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
        }
      }
    }
    break;
  }
  case BRACKET_ARAKAWA_OLD: {
    // Arakawa scheme for perpendicular flow
    
    result.allocate();

    const int ncz = mesh->LocalNz;
    const BoutReal partialFactor = 1.0/(12 * metric->dz);
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++){
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++){
        const BoutReal spacingFactor = partialFactor / metric->dx(jx, jy);
        for(int jz=0;jz<ncz;jz++) {
	  const int jzp = jz+1 < ncz ? jz + 1 : 0;
	  //Above is alternative to const int jzp = (jz + 1) % ncz;
	  const int jzm = jz-1 >=  0 ? jz - 1 : ncz-1;
	  //Above is alternative to const int jzm = (jz - 1 + ncz) % ncz;
          
          // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
          BoutReal Jpp = ( (f(jx,jy,jzp) - f(jx,jy,jzm))*
			   (g(jx+1,jy,jz) - g(jx-1,jy,jz)) -
			   (f(jx+1,jy,jz) - f(jx-1,jy,jz))*
			   (g(jx,jy,jzp) - g(jx,jy,jzm)) );

          // J+x
          BoutReal Jpx = ( g(jx+1,jy,jz)*(f(jx+1,jy,jzp)-f(jx+1,jy,jzm)) -
			   g(jx-1,jy,jz)*(f(jx-1,jy,jzp)-f(jx-1,jy,jzm)) -
			   g(jx,jy,jzp)*(f(jx+1,jy,jzp)-f(jx-1,jy,jzp)) +
			   g(jx,jy,jzm)*(f(jx+1,jy,jzm)-f(jx-1,jy,jzm)));

          // Jx+
          BoutReal Jxp = ( g(jx+1,jy,jzp)*(f(jx,jy,jzp)-f(jx+1,jy,jz)) -
			   g(jx-1,jy,jzm)*(f(jx-1,jy,jz)-f(jx,jy,jzm)) -
			   g(jx-1,jy,jzp)*(f(jx,jy,jzp)-f(jx-1,jy,jz)) +
			   g(jx+1,jy,jzm)*(f(jx+1,jy,jz)-f(jx,jy,jzm)));
          
          result(jx,jy,jz) = (Jpp + Jpx + Jxp) * spacingFactor;
        }
      }
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
    result = b0xGrad_dot_Grad(f, g) / metric->Bxy;
  }
  }
  
  result.setLocation(result_loc) ;
  
  return result;
}
