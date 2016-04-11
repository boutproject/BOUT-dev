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

const Field2D Grad_par(const Field2D &var, DIFF_METHOD method, CELL_LOC outloc) {
	return Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, CELL_LOC outloc, DIFF_METHOD method) {
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

const Field3D Grad_par(const Field3D &var, DIFF_METHOD method, CELL_LOC outloc) {
	return Grad_par(var, outloc, method);
}

// Model dvar/dt = Grad_par(f) with a maximum velocity of Vmax
const Field3D Grad_par(const Field3D &f, const Field3D &var, const Field2D &Vmax) {
#ifdef CHECK
	int msg_pos = msg_stack.push("Grad_par( Field3D, Field3D, Field2D )");
#endif
	Field2D sg = sqrt(mesh->g_22);
	Field3D result = DDY_MUSCL(f, var, sg*Vmax)/sg;
#ifdef CHECK
	msg_stack.pop(msg_pos);
#endif
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
  
	Field3D gys;
	gys.allocate();

	// Need Y derivative everywhere
	for(int x=1;x<=mesh->ngx-2;x++)
		for(int y=1;y<=mesh->ngy-2;y++)
	for(int z=0;z<ncz;z++) {
		gys(x, y, z) = (f(x, y+1, z) - f(x, y-1, z))/(0.5*mesh->dy(x, y+1) + mesh->dy(x, y) + 0.5*mesh->dy(x, y-1));
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
			BoutReal by = 1./sqrt(mesh->g_22(x, y));
			for(int z=0;z<ncz;z++) {
				int zm = (z - 1 + ncz) % ncz;
				int zp = (z + 1) % ncz;
        
				// bx = -DDZ(apar)
				bx(x, y, z) = (as(x, y, zm) - as(x, y, zp))/(2.*mesh->dz);
				// bz = DDX(f)
				bz(x, y, z) = (as(x+1, y, z) - as(x-1, y, z))/(0.5*mesh->dx(x-1, y) + mesh->dx(x, y) + 0.5*mesh->dx(x+1, y));
        
				// Now calculate (bx*d/dx + by*d/dy + bz*d/dz) f
        
				// Length dl for predictor
				BoutReal dl = fabs(mesh->dx(x, y)) / (fabs(bx(x, y, z)) + 1e-16);
				dl = BOUTMIN(dl, fabs(mesh->dy(x, y)) / (fabs(by) + 1e-16));
				dl = BOUTMIN(dl, mesh->dz / (fabs(bz(x, y, z)) + 1e-16));
        
				BoutReal fp, fm;
        
				// X differencing
				fp = fs(x+1, y, z)
					+ (0.25*dl/mesh->dz) * bz(x, y, z) * (fs(x+1, y, zm) - fs(x+1, y, zp))
						- 0.5*dl * by * gys(x+1, y, z);
        
				fm = fs(x-1, y, z)
					+ (0.25*dl/mesh->dz) * bz(x, y, z) * (fs(x-1, y, zm) - fs(x-1, y, zp))
						- 0.5*dl * by * gys(x-1, y, z);
        
				result(x, y, z) = bx(x, y, z) * (fp - fm) / (0.5*mesh->dx(x-1, y) + mesh->dx(x, y) + 0.5*mesh->dx(x+1, y));

				// Z differencing
        
				fp = fs(x, y, zp)
					+ (0.25*dl/mesh->dx(x, y)) * bx(x, y, z) * (fs(x-1, y, zp) - fs(x+1, y, zp))
						- 0.5*dl * by * gys(x, y, zp);
        
				fm = fs(x, y, zm)
					+ (0.25*dl/mesh->dx(x, y)) * bx(x, y, z) * (fs[x-1][y][zm] - fs(x+1, y, zm))
						- 0.5*dl * by * gys(x, y, zm);

				result(x, y, z) += bz(x, y, z) * (fp - fm) / (2.*mesh->dz);
        
				// Y differencing. Need X derivative
				gx(x, y, z) = (fs(x+1, y, z) - fs(x-1, y, z))/(0.5*mesh->dx(x-1, y) + mesh->dx(x, y) + 0.5*mesh->dx(x+1, y));
        
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
			BoutReal by = 1./sqrt(mesh->g_22(x, y));
			for(int z=0;z<ncz;z++) {
				int zm = (z - 1 + ncz) % ncz;
				int zp = (z + 1) % ncz;
        
				// Re-calculate dl
				BoutReal dl = fabs(mesh->dx(x, y)) / (fabs(bx(x, y, z)) + 1e-16);
				dl = BOUTMIN(dl, fabs(mesh->dy(x, y)) / (fabs(by) + 1e-16));
				dl = BOUTMIN(dl, mesh->dz / (fabs(bz(x, y, z)) + 1e-16));
        
				BoutReal fp, fm;

				fp = f[x][y+1][z]
					- 0.5*dl * bx[x][y][z] * gx[x][y+1][z]
						+ (0.25*dl/mesh->dz)   * bz[x][y][z] * (fs[x][y+1][zm] - fs[x][y+1][zp]);
        
				fm = f[x][y-1][z]
					- 0.5*dl * bx[x][y][z] * gx[x][y-1][z]
						+ (0.25*dl/mesh->dz)   * bz[x][y][z] * (fs[x][y-1][zm] - fs[x][y-1][zp]);

				result[x][y][z] += by * (fp - fm) / (0.5*mesh->dy[x][y-1] + mesh->dy[x][y] + 0.5*mesh->dy[x][y+1]);
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
	return VDDY(v, f)/sqrt(mesh->g_22);
}

const Field3D Vpar_Grad_par(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
	return VDDY(v, f, outloc, method)/sqrt(mesh->g_22);
}

const Field3D Vpar_Grad_par(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc) {
	return Vpar_Grad_par(v, f, outloc, method);
}

/*******************************************************************************
* Div_par
* parallel divergence operator B \partial_{||} (F/B)
*******************************************************************************/

const Field2D Div_par(const Field2D &f) {
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

const Field3D Div_par(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
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

const Field3D Div_par(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
	return Div_par(f, outloc, method);
}

//////// Flux methods

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
	return -mesh->Bxy*FDDY(v, f/mesh->Bxy, outloc, method)/sqrt(mesh->g_22);
}

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc) {
	return Div_par_flux(v,f, outloc, method);
}

//////// MUSCL schemes

const Field3D Div_par(const Field3D &f, const Field3D &var, const Field2D &Vmax) {
#ifdef CHECK
	int msg_pos = msg_stack.push("Div_par( Field3D, Field3D, Field2D )");
#endif

	Field3D result = mesh->Bxy*Grad_par(f/mesh->Bxy, var, Vmax/mesh->Bxy);
  
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
				result(jx, jy, jz) = 2.*(var(jx, jy, jz) - var(jx, jy-1, jz)) / (mesh->dy(jx, jy) * sqrt(mesh->g_22(jx, jy)) + mesh->dy(jx, jy-1) * sqrt(mesh->g_22(jx, jy-1)));
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
  
	for(int jx=0; jx<mesh->ngx;jx++) {
		for(int jy=0;jy<mesh->ngy-1;jy++) {
			for(int jz=0;jz<mesh->ngz;jz++) {
				result(jx, jy, jz) = (var(jx, jy+1, jz) - var(jx, jy, jz)) / (mesh->dy(jx, jy) * sqrt(mesh->g_22(jx, jy)));
			}
		}
	}
  
	/*
	bindex bx;
	bstencil f;
  
	result.allocate();

	start_index(&bx);
	do {
	var.setStencil(&f, &bx);
    
	result(bx.jx, bx.jy, bx.jz) = (f.yp - f.cc) / (mesh->dy(bx.jx, bx.jy) * sqrt(mesh->g_22(bx.jx, bx.jy)));
	}while(next_index3(&bx));
	*/
  
	return result;
}

const Field3D Div_par_LtoC(const Field2D &var) {
	return mesh->Bxy*Grad_par_LtoC(var/mesh->Bxy);
}

const Field3D Div_par_LtoC(const Field3D &var) {
	return mesh->Bxy*Grad_par_LtoC(var/mesh->Bxy);
}

const Field3D Div_par_CtoL(const Field2D &var) {
	return mesh->Bxy*Grad_par_CtoL(var/mesh->Bxy);
}

const Field3D Div_par_CtoL(const Field3D &var) {
	return mesh->Bxy*Grad_par_CtoL(var/mesh->Bxy);
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

const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc) {

#ifdef CHECK
	int msg_pos = msg_stack.push("Grad2_par2( Field3D )");
#endif

	Field2D sg;
	Field3D result, r2;
///#pragma omp parallel sections
///	{
///#pragma omp section
///		{
			sg = sqrt(mesh->g_22);
			sg = DDY(1./sg) / sg;
///		}
		if (sg.getLocation() != outloc) {
			mesh->communicate(sg);
			sg = interp_to(sg, outloc);
		}
///    
///#pragma omp section
		result = DDY(f,outloc);
    
///#pragma omp section
		r2 = D2DY2(f,outloc)/interp_to(mesh->g_22,outloc);
///	}
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
#ifdef CHECK
	int msg_pos = msg_stack.push("Delp2( Field2D )");
#endif

	Field2D result =  mesh->G1*DDX(f) + mesh->g11*D2DX2(f);

#ifdef CHECK
	msg_stack.pop(msg_pos);
#endif

	return result;
}

const Field3D Delp2(const Field3D &f, BoutReal zsmooth) {
	int msg_pos = msg_stack.push("Delp2( Field3D )");

	//return mesh->G1*DDX(f) + mesh->G3*DDZ(f) + mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f); //+ 2.0*mesh->g13*D2DXDZ(f)

	ASSERT2(mesh->xstart > 0); // Need at least one guard cell
  
	Field3D result;
	result.allocate();

	BoutReal ***fd, ***rd;
	fd = f.getData();
	rd = result.getData();

	int ncz = mesh->ngz-1;
  
	static dcomplex **ft = (dcomplex**) NULL, **delft;
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
			ZFFT(fd[jx][jy], mesh->zShift(jx, jy), ft[jx]);

		Laplacian::defaultInstance();
		// Loop over kz
#pragma omp parallel for
		for(int jz=0;jz<=ncz/2;jz++) {
			BoutReal filter;
			dcomplex a, b, c;
      
			if ((zsmooth > 0.0) && (jz > (int) (zsmooth*((BoutReal) ncz)))) filter=0.0; else filter=1.0;

			// No smoothing in the x direction
			for(int jx=mesh->xstart;jx<=mesh->xend;jx++) {
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
		for(int jx=mesh->xstart;jx<=mesh->xend;jx++) {

			ZFFT_rev(delft[jx], mesh->zShift[jx][jy], rd[jx][jy]);
			rd[jx][jy][ncz] = rd[jx][jy][0];
		}

		// Boundaries
		for(int jz=0;jz<ncz;jz++) {
			rd[0][jy][jz] = 0.0;
			rd[mesh->ngx-1][jy][jz] = 0.0;
		}
	}
  
	msg_stack.pop(msg_pos);

	// Set the output location
	result.setLocation(f.getLocation());

	return result;
}

const FieldPerp Delp2(const FieldPerp &f, BoutReal zsmooth) {
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
		ZFFT(fd[jx], mesh->zShift(jx, jy), ft[jx]);

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
	return D2DY2(f)/mesh->g_22 + DDY(mesh->J/mesh->g_22)*DDY(f)/mesh->J;
}

const Field3D Laplace_par(const Field3D &f) {
	return D2DY2(f)/mesh->g_22 + DDY(mesh->J/mesh->g_22)*DDY(f)/mesh->J;
}

/*******************************************************************************
* Laplacian
* Full Laplacian operator on scalar field
*******************************************************************************/

const Field2D Laplace(const Field2D &f) {
#ifdef CHECK
	int msg_pos = msg_stack.push("Laplace( Field2D )");
#endif

	Field2D result =  mesh->G1*DDX(f) + mesh->G2*DDY(f)
		+ mesh->g11*D2DX2(f) + mesh->g22*D2DY2(f)
			+ 2.0*mesh->g12*D2DXDY(f);

#ifdef CHECK
	msg_stack.pop(msg_pos);
#endif

	return result;
}

const Field3D Laplace(const Field3D &f) {
#ifdef CHECK
	int msg_pos = msg_stack.push("Laplace( Field3D )");
#endif

	Field3D result  = mesh->G1*DDX(f) + mesh->G2*DDY(f) + mesh->G3*DDZ(f)
		+ mesh->g11*D2DX2(f) + mesh->g22*D2DY2(f) + mesh->g33*D2DZ2(f)
			+ 2.0*(mesh->g12*D2DXDY(f) + mesh->g13*D2DXDZ(f) + mesh->g23*D2DYDZ(f));

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

const Field2D b0xGrad_dot_Grad(const Field2D &phi, const Field2D &A) {
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

const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A) {
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

const Field3D b0xGrad_dot_Grad(const Field3D &p, const Field2D &A, CELL_LOC outloc) {
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

const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc) {
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

const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method, CELL_LOC outloc, Solver *solver) {
	Field2D result;
	
	// Sort out cell locations
	CELL_LOC f_loc = f.getLocation(); // Location of f
	CELL_LOC g_loc = g.getLocation(); // Location of g
	CELL_LOC result_loc ; 
	  
	if(mesh->StaggerGrids) {
		if(outloc == CELL_DEFAULT){
			// Check that f and g are in the same location
			if (f_loc != g_loc){
				throw BoutException("Bracket currently requires both fields to have the same cell location");
			}
			else{
				result_loc = f_loc;      	  // Location of result
			}
		}
		else if (outloc != CELL_DEFAULT){
			// Check that f, and g are in the same location as the specified output location
			if(f_loc != g_loc || f_loc != outloc){
				throw BoutException("Bracket currently requires the location of both fields and the output locaton to be the same");
			}
			else{
				result_loc = outloc;      	  // Location of result
			}
		}
	}
	else{
		result_loc = CELL_CENTRE;      	  // Location of result		
	}
	
	if( (method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
		// Use a subset of terms for comparison to BOUT-06
		result = 0.0;
	}else {
		// Use full expression with all terms
		result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
	}
	
	result.setLocation(result_loc) ;
	
	return result;
}

const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method, CELL_LOC outloc, Solver *solver) {
	Field3D result;
	
	// Sort out cell locations
	CELL_LOC f_loc = f.getLocation(); // Location of f
	CELL_LOC g_loc = g.getLocation(); // Location of g
	CELL_LOC result_loc ; 
	  
	if(mesh->StaggerGrids) {
		if(outloc == CELL_DEFAULT){
			// Check that f and g are in the same location
			if (f_loc != g_loc){
				throw BoutException("Bracket currently requires both fields to have the same cell location");
			}
			else{
				result_loc = f_loc;      	  // Location of result
			}
		}
		else if (outloc != CELL_DEFAULT){
			// Check that f, and g are in the same location as the specified output location
			if(f_loc != g_loc || f_loc != outloc){
				throw BoutException("Bracket currently requires the location of both fields and the output locaton to be the same");
			}
			else{
				result_loc = outloc;      	  // Location of result
			}
		}
	}
	else{
		result_loc = CELL_CENTRE;      	  // Location of result		
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
    
			int ncz = mesh->ngz - 1;
			for(int x=mesh->xstart;x<=mesh->xend;x++)
			for(int y=mesh->ystart;y<=mesh->yend;y++) {
				for(int z=0;z<ncz;z++) {
					int zm = (z - 1 + ncz) % ncz;
					int zp = (z + 1) % ncz;
          
					BoutReal gp, gm;

					// Vx = DDZ(f)
					BoutReal vx = (f[x][y][zp] - f[x][y][zm])/(2.*mesh->dz);
          
					// Set stability condition
					solver->setMaxTimestep(mesh->dx[x][y] / (fabs(vx) + 1e-16));
          
					// X differencing
					if(vx > 0.0) {
						gp = g[x][y];
            
						gm = g[x-1][y];
            
					}else {
						gp = g[x+1][y];
            
						gm = g[x][y];
					}
          
					result[x][y][z] = vx * (gp - gm) / mesh->dx[x][y];
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
								/ (mesh->dx[jx][jy] * mesh->dz);

				// J+x
				BoutReal Jpx = 0.25*( g[jx+1][jy]*(fs[jx+1][jy][jzp]-fs[jx+1][jy][jzm]) -
					g[jx-1][jy]*(fs[jx-1][jy][jzp]-fs[jx-1][jy][jzm]) -
						g[jx][jy]*(fs[jx+1][jy][jzp]-fs[jx-1][jy][jzp]) +
							g[jx][jy]*(fs[jx+1][jy][jzm]-fs[jx-1][jy][jzm]))
								/ (mesh->dx[jx][jy] * mesh->dz);
				// Jx+
				BoutReal Jxp = 0.25*( g[jx+1][jy]*(fs[jx][jy][jzp]-fs[jx+1][jy][jz]) -
					g[jx-1][jy]*(fs[jx-1][jy][jz]-fs[jx][jy][jzm]) -
						g[jx-1][jy]*(fs[jx][jy][jzp]-fs[jx-1][jy][jz]) +
							g[jx+1][jy]*(fs[jx+1][jy][jz]-fs[jx][jy][jzm]))
								/ (mesh->dx[jx][jy] * mesh->dz);
          
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
			result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
		}
	}
	result.setLocation(result_loc) ;
	
	return result;
}

const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method, CELL_LOC outloc, Solver *solver) {
	Field3D result;
	
	// Sort out cell locations
	CELL_LOC f_loc = f.getLocation(); // Location of f
	CELL_LOC g_loc = g.getLocation(); // Location of g
	CELL_LOC result_loc ; 
	  
	if(mesh->StaggerGrids) {
		if(outloc == CELL_DEFAULT){
			// Take care of CELL_DEFAULT case
			// Check that f and g are in the same location
			if (f_loc != g_loc){
				throw BoutException("Bracket currently requires both fields to have the same cell location");
			}
			else{
				result_loc = f_loc;      	  // Location of result
			}
		}
		else if (outloc != CELL_DEFAULT){
			// Check that f, and g are in the same location as the specified output location
			if(f_loc != g_loc || f_loc != outloc){
				throw BoutException("Bracket currently requires the location of both fields and the output locaton to be the same");
			}
			else{
				result_loc = outloc;      	  // Location of result
			}
		}
	}
	else{
		result_loc = CELL_CENTRE;      	  // Location of result		
	}
	
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
			result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
		}
	}
	
	result.setLocation(result_loc) ;
	
	return result;
}

const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method, CELL_LOC outloc, Solver *solver) {
	
	Field3D result;
	
	// Sort out cell locations
	CELL_LOC f_loc = f.getLocation(); // Location of f
	CELL_LOC g_loc = g.getLocation(); // Location of g
	CELL_LOC result_loc ; 
	  
	if(mesh->StaggerGrids) {
		if(outloc == CELL_DEFAULT){
			// Take care of CELL_DEFAULT case
			// Check that f and g are in the same location
			if (f_loc != g_loc){
				throw BoutException("Bracket currently requires both fields to have the same cell location");
			}
			else{
				result_loc = f_loc;      	  // Location of result
			}
		}
		else if (outloc != CELL_DEFAULT){
			// Check that f, and g are in the same location as the specified output location
			if(f_loc != g_loc || f_loc != outloc){
				throw BoutException("Bracket currently requires the location of both fields and the output locaton to be the same");
			}
			else{
				result_loc = outloc;      	  // Location of result
			}
		}
	}
	else{
		result_loc = CELL_CENTRE;      	  // Location of result		
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
						vx[x][z] = (fs[x][y][zp] - fs[x][y][zm])/(2.*mesh->dz);
						// Vz = -DDX(f)
						vz[x][z] = (fs[x-1][y][z] - fs[x+1][y][z])/(0.5*mesh->dx[x-1][y] + mesh->dx[x][y] + 0.5*mesh->dx[x+1][y]);
          
						// Set stability condition
						solver->setMaxTimestep(fabs(mesh->dx[x][y]) / (fabs(vx[x][z]) + 1e-16));
						solver->setMaxTimestep(mesh->dz / (fabs(vz[x][z]) + 1e-16));
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
						gp = gs[x][y][z]
							+ (0.5*dt/mesh->dz) * ( (vz[x][z] > 0) ? vz[x][z]*(gs[x][y][zm] - gs[x][y][z]) : vz[x][z]*(gs[x][y][z] - gs[x][y][zp]) );
            
            
						gm = gs[x-1][y][z]
							//+ (0.5*dt/mesh->dz) * ( (vz[x-1][z] > 0) ? vz[x-1][z]*(g[x-1][y][zm] - g[x-1][y][z]) : vz[x-1][z]*(g[x-1][y][z] - g[x-1][y][zp]) );
							+ (0.5*dt/mesh->dz) * ( (vz[x][z] > 0) ? vz[x][z]*(gs[x-1][y][zm] - gs[x-1][y][z]) : vz[x][z]*(gs[x-1][y][z] - gs[x-1][y][zp]) );
            
					}else {
						gp = gs[x+1][y][z]
							//+ (0.5*dt/mesh->dz) * ( (vz[x+1][z] > 0) ? vz[x+1][z]*(gs[x+1][y][zm] - gs[x+1][y][z]) : vz[x+1][z]*(gs[x+1][y][z] - gs[x+1][y][zp]) );
							+ (0.5*dt/mesh->dz) * ( (vz[x][z] > 0) ? vz[x][z]*(gs[x+1][y][zm] - gs[x+1][y][z]) : vz[x][z]*(gs[x+1][y][z] - gs[x+1][y][zp]) );
            
						gm = gs[x][y][z] 
							+ (0.5*dt/mesh->dz) * ( (vz[x][z] > 0) ? vz[x][z]*(gs[x][y][zm] - gs[x][y][z]) : vz[x][z]*(gs[x][y][z] - gs[x][y][zp]) );
					}
          
					result[x][y][z] = vx[x][z] * (gp - gm) / mesh->dx[x][y];
          
					// Z differencing
					if(vz[x][z] > 0.0) {
						gp = gs[x][y][z]
							+ (0.5*dt/mesh->dx[x][y]) * ( (vx[x][z] > 0) ? vx[x][z]*(gs[x-1][y][z] - gs[x][y][z]) : vx[x][z]*(gs[x][y][z] - gs[x+1][y][z]) );
            
						gm = gs[x][y][zm]
							//+ (0.5*dt/mesh->dx[x][y]) * ( (vx[x][zm] > 0) ? vx[x][zm]*(gs[x-1][y][zm] - gs[x][y][zm]) : vx[x][zm]*(gs[x][y][zm] - gs[x+1][y][zm]) );
							+ (0.5*dt/mesh->dx[x][y]) * ( (vx[x][z] > 0) ? vx[x][z]*(gs[x-1][y][zm] - gs[x][y][zm]) : vx[x][z]*(gs[x][y][zm] - gs[x+1][y][zm]) );
					}else {
						gp = gs[x][y][zp]
							//+ (0.5*dt/mesh->dx[x][y]) * ( (vx[x][zp] > 0) ? vx[x][zp]*(gs[x-1][y][zp] - gs[x][y][zp]) : vx[x][zp]*(gs[x][y][zp] - gs[x+1][y][zp]) );
							+ (0.5*dt/mesh->dx[x][y]) * ( (vx[x][z] > 0) ? vx[x][z]*(gs[x-1][y][zp] - gs[x][y][zp]) : vx[x][z]*(gs[x][y][zp] - gs[x+1][y][zp]) );
            
						gm = gs[x][y][z]
							+ (0.5*dt/mesh->dx[x][y]) * ( (vx[x][z] > 0) ? vx[x][z]*(gs[x-1][y][z] - gs[x][y][z]) : vx[x][z]*(gs[x][y][z] - gs[x+1][y][z]) );
					}
          
					result[x][y][z] += vz[x][z] * (gp - gm) / mesh->dz;
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
								/ (mesh->dx[jx][jy] * mesh->dz);

				// J+x
				BoutReal Jpx = 0.25*( gs[jx+1][jy][jz]*(fs[jx+1][jy][jzp]-fs[jx+1][jy][jzm]) -
					gs[jx-1][jy][jz]*(fs[jx-1][jy][jzp]-fs[jx-1][jy][jzm]) -
						gs[jx][jy][jzp]*(fs[jx+1][jy][jzp]-fs[jx-1][jy][jzp]) +
							gs[jx][jy][jzm]*(fs[jx+1][jy][jzm]-fs[jx-1][jy][jzm]))
								/ (mesh->dx[jx][jy] * mesh->dz);
				// Jx+
				BoutReal Jxp = 0.25*( gs[jx+1][jy][jzp]*(fs[jx][jy][jzp]-fs[jx+1][jy][jz]) -
					gs[jx-1][jy][jzm]*(fs[jx-1][jy][jz]-fs[jx][jy][jzm]) -
						gs[jx-1][jy][jzp]*(fs[jx][jy][jzp]-fs[jx-1][jy][jz]) +
							gs[jx+1][jy][jzm]*(fs[jx+1][jy][jz]-fs[jx][jy][jzm]))
								/ (mesh->dx[jx][jy] * mesh->dz);
          
				result[jx][jy][jz] = (Jpp + Jpx + Jxp) / 3.;
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
			result = b0xGrad_dot_Grad(f, g) / mesh->Bxy;
		}
	}
	
	result.setLocation(result_loc) ;
	
	return result;
}
