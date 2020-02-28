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
#include <bout/solver.hxx>
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

#include <cmath>

/*******************************************************************************
* Grad_par
* The parallel derivative along unperturbed B-field
*******************************************************************************/

const Coordinates::metric_field_type Grad_par(const Field2D& var, CELL_LOC outloc,
                                              const std::string& method) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
}

const Coordinates::metric_field_type
Grad_par(const Field2D& var, const std::string& method, CELL_LOC outloc) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, CELL_LOC outloc, const std::string &method) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, const std::string &method, CELL_LOC outloc) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
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
  ASSERT1_FIELDS_COMPATIBLE(apar, f);
  ASSERT1(f.hasParallelSlices());

  Mesh *mesh = apar.getMesh();

  Field3D result{emptyFrom(f)};
  
  int ncz = mesh->LocalNz;

  Coordinates *metric = apar.getCoordinates();

  Field3D gys{emptyFrom(f)};

  // Need Y derivative everywhere
  for(int x=1;x<=mesh->LocalNx-2;x++)
    for(int y=1;y<=mesh->LocalNy-2;y++)
      for(int z=0;z<ncz;z++) {
        gys(x, y, z) = (f.yup()(x, y+1, z) - f.ydown()(x, y-1, z))/(0.5*metric->dy(x, y+1, z) + metric->dy(x, y, z) + 0.5*metric->dy(x, y-1, z));
      }
  
  for(int x=1;x<=mesh->LocalNx-2;x++) {
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      for(int z=0;z<ncz;z++) {
	BoutReal by = 1./sqrt(metric->g_22(x, y, z));
        // Z indices zm and zp
        int zm = (z - 1 + ncz) % ncz;
        int zp = (z + 1) % ncz;
        
        // bx = -DDZ(apar)
        BoutReal bx = (apar(x, y, zm) - apar(x, y, zp))/(0.5*metric->dz(x, y, zm) + metric->dz(x, y, z) + 0.5*metric->dz(x, y, zp));
        // bz = DDX(f)
        BoutReal bz = (apar(x+1, y, z) - apar(x-1, y, z))/(0.5*metric->dx(x-1, y, z) + metric->dx(x, y, z) + 0.5*metric->dx(x+1, y, z));
        
	// Now calculate (bx*d/dx + by*d/dy + bz*d/dz) f
        
        // Length dl for predictor
        BoutReal dl = fabs(metric->dx(x, y, z)) / (fabs(bx) + 1e-16);
        dl = BOUTMIN(dl, fabs(metric->dy(x, y, z)) / (fabs(by) + 1e-16));
        dl = BOUTMIN(dl, fabs(metric->dz(x, y, z)) / (fabs(bz) + 1e-16));
        
	BoutReal fp, fm;
        
        // X differencing
        fp = f(x+1, y, z)
          + (0.25*dl/metric->dz(x, y, z)) * bz * (f(x+1, y, zm) - f(x+1, y, zp))
          - 0.5*dl * by * gys(x+1, y, z);
        
        fm = f(x-1, y, z)
          + (0.25*dl/metric->dz(x, y, z)) * bz * (f(x-1, y, zm) - f(x-1, y, zp))
          - 0.5*dl * by * gys(x-1, y, z);
        
        result(x, y, z) = bx * (fp - fm) / (0.5*metric->dx(x-1, y, z) + metric->dx(x, y, z) + 0.5*metric->dx(x+1, y, z));

	// Z differencing
        
        fp = f(x, y, zp)
          + (0.25*dl/metric->dx(x, y, z)) * bx * (f(x-1, y, zp) - f(x+1, y, zp))
          - 0.5*dl * by * gys(x, y, zp);
        
        fm = f(x, y, zm)
          + (0.25*dl/metric->dx(x, y, z)) * bx * (f(x-1,y,zm) - f(x+1, y, zm))
          - 0.5*dl * by * gys(x, y, zm);

        result(x, y, z) += bz * (fp - fm) / (0.5*metric->dz(x, y, zm) + metric->dz(x, y, z) + 0.5*metric->dz(x, y, zp));

        // Y differencing
        
        fp = f.yup()(x,y+1,z)
          - 0.5*dl * bx * (f.yup()(x+1, y+1, z) - f.yup()(x-1, y+1, z))/(0.5*metric->dx(x-1, y, z) + metric->dx(x, y, z) + 0.5*metric->dx(x+1, y, z))
          
          + (0.25*dl/metric->dz(x,y,z)) * bz * (f.yup()(x,y+1,zm) - f.yup()(x,y+1,zp));
        
        fm = f.ydown()(x,y-1,z)
          - 0.5*dl * bx * (f.ydown()(x+1, y-1, z) - f.ydown()(x-1, y-1, z))/(0.5*metric->dx(x-1, y, z) + metric->dx(x, y, z) + 0.5*metric->dx(x+1, y, z))
          + (0.25*dl/metric->dz(x,y,z)) * bz * (f.ydown()(x,y-1,zm) - f.ydown()(x,y-1,zp));

        result(x,y,z) += by * (fp - fm) / (0.5*metric->dy(x,y-1, z) + metric->dy(x,y,z) + 0.5*metric->dy(x,y+1,z));
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

const Coordinates::metric_field_type Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                                   CELL_LOC outloc,
                                                   const std::string& method) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

const Coordinates::metric_field_type Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                                   const std::string& method,
                                                   CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, const std::string &method, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

/*******************************************************************************
* Div_par
* parallel divergence operator B \partial_{||} (F/B)
*******************************************************************************/
const Coordinates::metric_field_type Div_par(const Field2D& f, CELL_LOC outloc,
                                             const std::string& method) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Coordinates::metric_field_type Div_par(const Field2D& f, const std::string& method,
                                             CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D &f, const std::string &method, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D& f, const Field3D& v) {
  ASSERT1_FIELDS_COMPATIBLE(f, v);
  ASSERT1(f.hasParallelSlices());
  ASSERT1(v.hasParallelSlices());

  // Parallel divergence, using velocities at cell boundaries
  // Note: Not guaranteed to be flux conservative
  Mesh* mesh = f.getMesh();

  Field3D result{emptyFrom(f)};

  Coordinates* coord = f.getCoordinates();

  for (int i = mesh->xstart; i <= mesh->xend; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = mesh->zstart; k <= mesh->zend; k++) {
        // Value of f and v at left cell face
        BoutReal fL = 0.5 * (f(i, j, k) + f.ydown()(i, j - 1, k));
        BoutReal vL = 0.5 * (v(i, j, k) + v.ydown()(i, j - 1, k));

        BoutReal fR = 0.5 * (f(i, j, k) + f.yup()(i, j + 1, k));
        BoutReal vR = 0.5 * (v(i, j, k) + v.yup()(i, j + 1, k));

        // Calculate flux at right boundary (y+1/2)
        BoutReal fluxRight =
            fR * vR * (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        // Calculate at left boundary (y-1/2)
        BoutReal fluxLeft =
            fL * vL * (coord->J(i, j, k) + coord->J(i, j - 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        result(i, j, k) =
            (fluxRight - fluxLeft) / (coord->dy(i, j, k) * coord->J(i, j, k));
      }
    }

  return result;
}

//////// Flux methods

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method) {
  Coordinates *metric = f.getCoordinates(outloc);

  auto Bxy_floc = f.getCoordinates()->Bxy;

  if (!f.hasParallelSlices()) {
    return metric->Bxy*FDDY(v, f/Bxy_floc, outloc, method)/sqrt(metric->g_22);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  // Distinct fields
  f_B.splitParallelSlices();
  f_B.yup() = f.yup() / Bxy_floc;
  f_B.ydown() = f.ydown() / Bxy_floc;
  return metric->Bxy*FDDY(v, f_B, outloc, method)/sqrt(metric->g_22);
}

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, const std::string &method, CELL_LOC outloc) {
  return Div_par_flux(v,f, outloc, method);
}

/*******************************************************************************
* Grad2_par2
* second parallel derivative
*
* (b dot Grad)(b dot Grad)
*
* Note: For parallel Laplacian use LaplacePar
*******************************************************************************/

const Coordinates::metric_field_type Grad2_par2(const Field2D& f, CELL_LOC outloc,
                                                const std::string& method) {
  return f.getCoordinates(outloc)->Grad2_par2(f, outloc, method);
}

const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Grad2_par2(f, outloc, method);
}

/*******************************************************************************
* Div_par_K_Grad_par
* Parallel divergence of diffusive flux, K*Grad_par
*******************************************************************************/

const Coordinates::metric_field_type Div_par_K_Grad_par(BoutReal kY, const Field2D& f,
                                                        CELL_LOC outloc) {
  return kY*Grad2_par2(f, outloc);
}

const Field3D Div_par_K_Grad_par(BoutReal kY, const Field3D &f, CELL_LOC outloc) {
  return kY*Grad2_par2(f, outloc);
}

const Coordinates::metric_field_type
Div_par_K_Grad_par(const Field2D& kY, const Field2D& f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

const Field3D Div_par_K_Grad_par(const Field2D &kY, const Field3D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field2D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field3D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

/*******************************************************************************
* Delp2
* perpendicular Laplacian operator
*******************************************************************************/

const Coordinates::metric_field_type Delp2(const Field2D& f, CELL_LOC outloc, bool useFFT) {
  return f.getCoordinates(outloc)->Delp2(f, outloc, useFFT);
}

const Field3D Delp2(const Field3D& f, CELL_LOC outloc, bool useFFT) {
  return f.getCoordinates(outloc)->Delp2(f, outloc, useFFT);
}

const FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
  return f.getCoordinates(outloc)->Delp2(f, outloc, useFFT);
}

/*******************************************************************************
* LaplacePerp
* Full perpendicular Laplacian operator on scalar field
*
* Laplace_perp = Laplace - Laplace_par
*******************************************************************************/

const Coordinates::metric_field_type Laplace_perp(const Field2D& f, CELL_LOC outloc) {
  return Laplace(f, outloc) - Laplace_par(f, outloc);
}

const Field3D Laplace_perp(const Field3D &f, CELL_LOC outloc) {
  return Laplace(f, outloc) - Laplace_par(f, outloc);
}

/*******************************************************************************
* LaplacePar
* Full parallel Laplacian operator on scalar field
*
* LaplacePar(f) = Div( b (b dot Grad(f)) ) 
*
*******************************************************************************/

const Coordinates::metric_field_type Laplace_par(const Field2D& f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace_par(f, outloc);
}

const Field3D Laplace_par(const Field3D &f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace_par(f, outloc);
}

/*******************************************************************************
* Laplacian
* Full Laplacian operator on scalar field
*******************************************************************************/

const Coordinates::metric_field_type Laplace(const Field2D& f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace(f, outloc);
}

const Field3D Laplace(const Field3D &f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace(f, outloc);
}

/*******************************************************************************
* b0xGrad_dot_Grad
* Terms of form b0 x Grad(phi) dot Grad(A)
* Used for ExB terms and perturbed B field using A_||
*******************************************************************************/

const Coordinates::metric_field_type b0xGrad_dot_Grad(const Field2D& phi,
                                                      const Field2D& A, CELL_LOC outloc) {

  TRACE("b0xGrad_dot_Grad( Field2D , Field2D )");
  
  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(phi.getMesh() == A.getMesh());

  Coordinates *metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  Coordinates::metric_field_type dpdx = DDX(phi, outloc);
  Coordinates::metric_field_type dpdy = DDY(phi, outloc);

  // Calculate advection velocity
  Coordinates::metric_field_type vx = -metric->g_23 * dpdy;
  Coordinates::metric_field_type vy = metric->g_23 * dpdx;

  // Upwind A using these velocities
  Coordinates::metric_field_type result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);
  result /= metric->J*sqrt(metric->g_22);

  ASSERT1(result.getLocation() == outloc);

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field2D , Field3D )");

  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(phi.getMesh() == A.getMesh());

  Mesh *mesh = phi.getMesh();

  Coordinates *metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  Coordinates::metric_field_type dpdx = DDX(phi, outloc);
  Coordinates::metric_field_type dpdy = DDY(phi, outloc);

  // Calculate advection velocity
  Coordinates::metric_field_type vx = -metric->g_23 * dpdy;
  Coordinates::metric_field_type vy = metric->g_23 * dpdx;
  Coordinates::metric_field_type vz = metric->g_12 * dpdy - metric->g_22 * dpdx;

  if(mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion * vx;
  }

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc) + VDDZ(vz, A, outloc);

  result /= (metric->J*sqrt(metric->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif

  ASSERT2(result.getLocation() == outloc);
  
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &p, const Field2D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field3D , Field2D )");

  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(p.getMesh() == A.getMesh());

  Coordinates *metric = p.getCoordinates(outloc);

  // Calculate phi derivatives
  Field3D dpdx = DDX(p, outloc);
  Field3D dpdy = DDY(p, outloc);
  Field3D dpdz = DDZ(p, outloc);

  // Calculate advection velocity
  Field3D vx = metric->g_22 * dpdz - metric->g_23 * dpdy;
  Field3D vy = metric->g_23 * dpdx - metric->g_12 * dpdz;

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);

  result /= (metric->J*sqrt(metric->g_22));
  
#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+p.name+","+A.name+")";
#endif
  
  ASSERT2(result.getLocation() == outloc);

  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field3D , Field3D )");

  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(phi.getMesh() == A.getMesh());

  Mesh * mesh = phi.getMesh();

  Coordinates *metric = phi.getCoordinates(outloc);

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

  result /=  (metric->J*sqrt(metric->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif

  ASSERT2(result.getLocation() == outloc);

  return result;
}


const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, CELL_LOC outloc) {
  
  Field3D result = 0.0;

  //////////////////////////////////////////
  // X-Z diffusion
  // 
  //            Z
  //            |
  // 
  //     o --- gU --- o
  //     |     nU     |
  //     |            |
  //    gL nL      nR gR    -> X
  //     |            |
  //     |     nD     |
  //     o --- gD --- o
  //
  Coordinates *coords = a.getCoordinates(outloc);
  Mesh *mesh = f.getMesh();

  Field3D fs = f;
  Field3D as = a;

  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++) {

	// wrap k-index around as Z is (currently) periodic.
	int kp = (k+1) % (mesh->LocalNz);
	int km = (k-1+mesh->LocalNz) % (mesh->LocalNz);
        
        // Calculate gradients on cell faces -- assumes constant grid spacing
        
        BoutReal gR = (coords->g11(i,j,k) + coords->g11(i+1,j,k)) * (fs(i+1,j,k) - fs(i,j,k))/(coords->dx(i+1,j,k) + coords->dx(i,j,k))
          + 0.5*(coords->g13(i,j,k) + coords->g13(i+1,j,k))*(fs(i+1,j,kp) - fs(i+1,j,km) + fs(i,j,kp) - fs(i,j,km))/(4.*coords->dz(i,j,k));
        
	BoutReal gL = (coords->g11(i-1,j,k) + coords->g11(i,j,k))*(fs(i,j,k) - fs(i-1,j,k))/(coords->dx(i-1,j,k) + coords->dx(i,j,k))
	  + 0.5*(coords->g13(i-1,j,k) + coords->g13(i,j,k))*(fs(i-1,j,kp) - fs(i-1,j,km) + f(i,j,kp) - f(i,j,km))/(4*coords->dz(i,j,k));
	
	BoutReal gD = coords->g13(i,j,k)*(fs(i+1,j,km) - fs(i-1,j,km) + fs(i+1,j,k) - fs(i-1,j,k))/(4.*coords->dx(i,j,k))
	  + coords->g33(i,j,k)*(fs(i,j,k) - fs(i,j,km))/coords->dz(i,j,k);
        
        BoutReal gU = coords->g13(i,j,k)*(fs(i+1,j,kp) - fs(i-1,j,kp) + fs(i+1,j,k) - fs(i-1,j,k))/(4.*coords->dx(i,j,k))
          + coords->g33(i,j,k)*(fs(i,j,kp) - fs(i,j,k))/coords->dz(i,j,k);
          
        
        // Flow right
        BoutReal flux = gR * 0.25*(coords->J(i+1,j,k) + coords->J(i,j,k)) *(as(i+1,j,k) + as(i,j,k));
        result(i,j,k)   += flux / (coords->dx(i,j,k)*coords->J(i,j,k));
        
        // Flow left
        flux = gL * 0.25*(coords->J(i-1,j,k) + coords->J(i,j,k)) *(as(i-1,j,k) + as(i,j,k));
        result(i,j,k)   -= flux / (coords->dx(i,j,k)*coords->J(i,j,k));
        
        
        // Flow up
        flux = gU * 0.25*(coords->J(i,j,k) + coords->J(i,j,kp)) *(as(i,j,k) + as(i,j,kp));
        result(i,j,k) += flux / (coords->dz(i,j,k) * coords->J(i,j,k));

	// Flow down
        flux = gD * 0.25*(coords->J(i,j,km) + coords->J(i,j,k)) *(as(i,j,km) + as(i,j,k));
        result(i,j,k) += flux / (coords->dz(i,j,k) * coords->J(i,j,k));
      }

  return result;
}


/*******************************************************************************
 * Poisson bracket
 * Terms of form b0 x Grad(f) dot Grad(g) / B = [f, g]
 *******************************************************************************/

const Coordinates::metric_field_type bracket(const Field2D& f, const Field2D& g,
                                             BRACKET_METHOD method, CELL_LOC outloc,
                                             Solver* UNUSED(solver)) {
  TRACE("bracket(Field2D, Field2D)");

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Coordinates::metric_field_type result{emptyFrom(f)};

  if( (method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
    result.setLocation(outloc);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / f.getCoordinates(outloc)->Bxy;
  }
  return result;
}

const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("bracket(Field3D, Field2D)");

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Mesh *mesh = f.getMesh();

  Field3D result{emptyFrom(f).setLocation(outloc)};

  Coordinates *metric = f.getCoordinates(outloc);

  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");

#ifndef COORDINATES_USE_3D
    const int ncz = mesh->LocalNz;

    for(int x=mesh->xstart;x<=mesh->xend;x++)
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        for (int z = 0; z < ncz; z++) {
          int zm = (z - 1 + ncz) % ncz;
	  int zp = (z + 1) % ncz;
          
	  BoutReal gp, gm;
	  
          // Vx = DDZ(f)
          BoutReal vx = (f(x,y,zp) - f(x,y,zm))/(2.*metric->dz(x,y,z));
          
          // Set stability condition
          solver->setMaxTimestep(metric->dx(x, y, z) / (fabs(vx) + 1e-16));

          // X differencing
          if(vx > 0.0) {
            gp = g(x,y);
            
            gm = g(x-1,y);
            
          }else {
            gp = g(x+1,y);
            
            gm = g(x,y);
          }

          result(x, y, z) = vx * (gp - gm) / metric->dx(x, y, z);
        }
      }
#else
    throw BoutException("BRACKET_CTU not valid with 3D metrics yet.");
#endif
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow. Here as a test

#ifndef COORDINATES_USE_3D
    const int ncz = mesh->LocalNz;

    BOUT_FOR(j2D, result.getRegion2D("RGN_NOBNDRY")) {
      // Get constants for this iteration
      const BoutReal fac = 1.0 / (12 * metric->dz[j2D]);
      const BoutReal spacingFactor = fac / metric->dx[j2D];
      const int jy = j2D.y(), jx = j2D.x();
      const int xm = jx - 1, xp = jx + 1;

      // Extract relevant Field2D values
      const BoutReal gxm = g(xm, jy), gc = g(jx, jy), gxp = g(xp, jy);

      // Index Field3D as 2D to get start of z data block
      const auto fxm = f(xm, jy), fc = f(jx, jy), fxp = f(xp, jy);

      // Here we split the loop over z into three parts; the first value, the middle block
      // and the last value
      // this is to allow the loop used in the middle block to vectorise.

      // The first value
      {
        const int jzp = 1;
        const int jzm = ncz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = gxp * (fxp[jzp] - fxp[jzm]) - gxm * (fxm[jzp] - fxm[jzm]) +
                             gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]);

        result(jx, jy, 0) = (Jpp + Jpx) * spacingFactor;
      }

      // The middle block
      for (int jz = 1; jz < mesh->LocalNz - 1; jz++) {
        const int jzp = jz + 1;
        const int jzm = jz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = gxp * (fxp[jzp] - fxp[jzm]) - gxm * (fxm[jzp] - fxm[jzm]) +
                             gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]);

        result(jx, jy, jz) = (Jpp + Jpx) * spacingFactor;
      }

      // The last value
      {
        const int jzp = 0;
        const int jzm = ncz - 2;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = gxp * (fxp[jzp] - fxp[jzm]) - gxm * (fxm[jzp] - fxm[jzm]) +
                             gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]);

        result(jx, jy, ncz - 1) = (Jpp + Jpx) * spacingFactor;
      }
    }
#else
    throw BoutException("BRACKET_ARAKAWA not valid with 3D metrics yet.");
#endif

    break;
  }
  case BRACKET_ARAKAWA_OLD: {
#ifndef COORDINATES_USE_3D
    const int ncz = mesh->LocalNz;
    BOUT_OMP(parallel for)
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++){
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++){
	const BoutReal partialFactor = 1.0/(12 * metric->dz(jx,jy));
    	const BoutReal spacingFactor = partialFactor / metric->dx(jx,jy);
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
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
#else
    throw BoutException("BRACKET_ARAKAWA_OLD not valid with 3D metrics yet.");
#endif
      break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy;
  }
  }
  return result;
}

const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("bracket(Field2D, Field3D)");

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation())

  Mesh *mesh = f.getMesh();

  Field3D result(mesh);

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
    result = VDDZ(-DDX(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    Coordinates *metric = f.getCoordinates(outloc);
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy;
  }
  }
  
  return result;
}

const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("Field3D, Field3D");

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Mesh *mesh = f.getMesh();

  Field3D result{emptyFrom(f).setLocation(outloc)};

  Coordinates *metric = f.getCoordinates(outloc);

  if (mesh->GlobalNx == 1 || mesh->GlobalNz == 1) {
    result=0;
    result.setLocation(outloc);
    return result;
  }
  
  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
#ifndef COORDINATES_USE_3D
    if(!solver)
      throw BoutException("CTU method requires access to the solver");

    // Get current timestep
    BoutReal dt = solver->getCurrentTimestep();
    
    FieldPerp vx(mesh), vz(mesh);
    vx.allocate();
    vx.setLocation(outloc);
    vz.allocate();
    vz.setLocation(outloc);
    
    int ncz = mesh->LocalNz;
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      for(int x=1;x<=mesh->LocalNx-2;x++) {
        for (int z = 0; z < mesh->LocalNz; z++) {
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;
          
          // Vx = DDZ(f)
          vx(x,z) = (f(x,y,zp) - f(x,y,zm))/(2.*metric->dz(x,y,z));
          // Vz = -DDX(f)
          vz(x,z) = (f(x-1,y,z) - f(x+1,y,z))/(0.5*metric->dx(x-1,y) + metric->dx(x,y) + 0.5*metric->dx(x+1,y));
          
          // Set stability condition
          solver->setMaxTimestep(fabs(metric->dx(x,y)) / (fabs(vx(x,z)) + 1e-16));
          solver->setMaxTimestep(metric->dz(x,y) / (fabs(vz(x,z)) + 1e-16));
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
                 (0.5 * dt / metric->dz(x,y)) * ((vz(x, z) > 0)
                                                ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                                : vz(x, z) * (g(x, y, z) - g(x, y, zp)));

            gm = g(x - 1, y, z) +
                 (0.5 * dt / metric->dz(x,y)) *
                     ((vz(x, z) > 0) ? vz(x, z) * (g(x - 1, y, zm) - g(x - 1, y, z))
                                     : vz(x, z) * (g(x - 1, y, z) - g(x - 1, y, zp)));

          } else {
            gp = g(x + 1, y, z) +
                 (0.5 * dt / metric->dz(x,y)) *
                     ((vz(x, z) > 0) ? vz(x, z) * (g(x + 1, y, zm) - g(x + 1, y, z))
                                     : vz[x][z] * (g(x + 1, y, z) - g(x + 1, y, zp)));

            gm = g(x, y, z) +
                 (0.5 * dt / metric->dz(x,y)) * ((vz(x, z) > 0)
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

          result(x, y, z) += vz(x, z) * (gp - gm) / metric->dz(x,y);
        }
    }
#else
    throw BoutException("BRACKET_CTU not valid with 3D metrics yet.");
#endif
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow
#ifndef COORDINATES_USE_3D
    const int ncz = mesh->LocalNz;

    // We need to discard const qualifier in order to manipulate
    // storage array directly
    Field3D f_temp = f;
    Field3D g_temp = g;

    BOUT_FOR(j2D, result.getRegion2D("RGN_NOBNDRY")) {
      const BoutReal partialFactor = 1.0/(12 * metric->dz[j2D]);
      const BoutReal spacingFactor = partialFactor / metric->dx[j2D];
      const int jy = j2D.y(), jx = j2D.x();
      const int xm = jx - 1, xp = jx + 1;

      const auto Fxm = f_temp(xm, jy), Fx = f_temp(jx, jy), Fxp = f_temp(xp, jy);
      const auto Gxm = g_temp(xm, jy), Gx = g_temp(jx, jy), Gxp = g_temp(xp, jy);

      // Here we split the loop over z into three parts; the first value, the middle block
      // and the last value
      // this is to allow the loop used in the middle block to vectorise.

      {
        const int jz = 0;
        const int jzp = 1;
        const int jzm = ncz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = ((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]) -
                              (Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm]));

        // J+x
        const BoutReal Jpx =
            (Gxp[jz] * (Fxp[jzp] - Fxp[jzm]) - Gxm[jz] * (Fxm[jzp] - Fxm[jzm]) -
             Gx[jzp] * (Fxp[jzp] - Fxm[jzp]) + Gx[jzm] * (Fxp[jzm] - Fxm[jzm]));

        // Jx+
        const BoutReal Jxp =
            (Gxp[jzp] * (Fx[jzp] - Fxp[jz]) - Gxm[jzm] * (Fxm[jz] - Fx[jzm]) -
             Gxm[jzp] * (Fx[jzp] - Fxm[jz]) + Gxp[jzm] * (Fxp[jz] - Fx[jzm]));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }

      for (int jz = 1; jz < mesh->LocalNz - 1; jz++) {
        const int jzp = jz + 1;
        const int jzm = jz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = ((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]) -
                              (Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm]));

        // J+x
        const BoutReal Jpx =
            (Gxp[jz] * (Fxp[jzp] - Fxp[jzm]) - Gxm[jz] * (Fxm[jzp] - Fxm[jzm]) -
             Gx[jzp] * (Fxp[jzp] - Fxm[jzp]) + Gx[jzm] * (Fxp[jzm] - Fxm[jzm]));

        // Jx+
        const BoutReal Jxp =
            (Gxp[jzp] * (Fx[jzp] - Fxp[jz]) - Gxm[jzm] * (Fxm[jz] - Fx[jzm]) -
             Gxm[jzp] * (Fx[jzp] - Fxm[jz]) + Gxp[jzm] * (Fxp[jz] - Fx[jzm]));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }

      {
        const int jz = ncz - 1;
        const int jzp = 0;
        const int jzm = ncz - 2;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = ((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]) -
                              (Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm]));

        // J+x
        const BoutReal Jpx =
            (Gxp[jz] * (Fxp[jzp] - Fxp[jzm]) - Gxm[jz] * (Fxm[jzp] - Fxm[jzm]) -
             Gx[jzp] * (Fxp[jzp] - Fxm[jzp]) + Gx[jzm] * (Fxp[jzm] - Fxm[jzm]));

        // Jx+
        const BoutReal Jxp =
            (Gxp[jzp] * (Fx[jzp] - Fxp[jz]) - Gxm[jzm] * (Fxm[jz] - Fx[jzm]) -
             Gxm[jzp] * (Fx[jzp] - Fxm[jz]) + Gxp[jzm] * (Fxp[jz] - Fx[jzm]));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }
    }
#else
    throw BoutException("BRACKET_ARAKAWA not valid with 3D metrics yet.");
#endif

    break;
  }
  case BRACKET_ARAKAWA_OLD: {
    // Arakawa scheme for perpendicular flow
#ifndef COORDINATES_USE_3D

    const int ncz = mesh->LocalNz;

    // We need to discard const qualifier in order to manipulate
    // storage array directly
    Field3D f_temp = f;
    Field3D g_temp = g;

    BOUT_OMP(parallel for)
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++){
      for(int jy=mesh->ystart;jy<=mesh->yend;jy++){
	const BoutReal spacingFactor = 1.0 / (12 * metric->dz(jx,jy)
					      * metric->dx(jx, jy));
        const BoutReal *Fxm = f_temp(jx-1, jy);
        const BoutReal *Fx  = f_temp(jx,   jy);
        const BoutReal *Fxp = f_temp(jx+1, jy);
        const BoutReal *Gxm = g_temp(jx-1, jy);
        const BoutReal *Gx  = g_temp(jx,   jy);
        const BoutReal *Gxp = g_temp(jx+1, jy);
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
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
#else
    throw BoutException("BRACKET_ARAKAWA_OLD not valid with 3D metrics yet.");
#endif
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f, outloc), g, outloc) + VDDZ(-DDX(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy;
  }
  }
  
  return result;
}
