
/*******************************************************************************
 * Differential operators
 *
 * Changelog:
 *
 * 2009-01 Ben Dudson <bd512@york.ac.uk> 
 *    * Added two optional parameters which can be put in any order
 *      These determine the method to use (DIFF_METHOD)
 *      and CELL_LOC location of the result.
 *      Both of these options are defined in bout_types.hxx
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
 *******************************************************************************/

#ifndef __DIFOPS_H__
#define __DIFOPS_H__

#include "field3d.hxx"
#include "field2d.hxx"

#include "bout_types.hxx"

#include "bout/solver.hxx"

// Parallel derivative (central differencing)
const Field2D Grad_par(const Field2D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field2D Grad_par(const Field2D &var, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);

const Field3D Grad_par(const Field3D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field3D Grad_par(const Field3D &var, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);

// MUSCL schemes. Model dvar/dt = Grad_par(f) with a maximum velocity of Vmax
const Field3D Grad_par(const Field3D &f, const Field3D &var, const Field2D &Vmax);
const Field3D Grad_par(const Field3D &f, const Field3D &var, BoutReal Vmax);

// b0 dot Grad  -  (1/B)b0 x Grad(apar) dot Grad
const Field3D Grad_parP(const Field3D &apar, const Field3D &f);

// vpar times parallel derivative (upwinding)
const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f);
const Field3D Vpar_Grad_par(const Field &v, const Field &f, 
			    CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field3D Vpar_Grad_par(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);


// parallel divergence operator B \partial_{||} (F/B)
const Field2D Div_par(const Field2D &f);
const Field3D Div_par(const Field3D &f, 
		      CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field3D Div_par(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT);

// Flux methods. Model divergence of flux: df/dt =  Div(v * f)
const Field3D Div_par_flux(const Field3D &v, const Field3D &f, 
		      CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field3D Div_par_flux(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT);

// MUSCL scheme. Model dvar/dt = Div_par(f) with a maximum velocity of Vmax
const Field3D Div_par(const Field3D &f, const Field3D &var, const Field2D &Vmax);
const Field3D Div_par(const Field3D &f, const Field3D &var, BoutReal Vmax);


// second parallel derivative
const Field2D Grad2_par2(const Field2D &f);
const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);

// Parallel derivatives, converting between cell-centred and lower cell boundary
// These are a simple way to do staggered differencing
const Field3D Grad_par_CtoL(const Field3D &var);
const Field3D Vpar_Grad_par_LCtoC(const Field &v, const Field &f);
const Field3D Grad_par_LtoC(const Field3D &var);
const Field3D Div_par_LtoC(const Field2D &var);
const Field3D Div_par_LtoC(const Field3D &var);
const Field3D Div_par_CtoL(const Field2D &var);
const Field3D Div_par_CtoL(const Field3D &var);

// Parallel divergence of diffusive flux, K*Grad_par
const Field2D Div_par_K_Grad_par(BoutReal kY, Field2D &f);
const Field3D Div_par_K_Grad_par(BoutReal kY, Field3D &f);
const Field2D Div_par_K_Grad_par(Field2D &kY, Field2D &f);
const Field3D Div_par_K_Grad_par(Field2D &kY, Field3D &f);
const Field3D Div_par_K_Grad_par(Field3D &kY, Field2D &f);
const Field3D Div_par_K_Grad_par(Field3D &kY, Field3D &f);

// Divergence of perpendicular diffusive flux kperp*Grad_perp
const Field3D Div_K_perp_Grad_perp(const Field2D &kperp, const Field3D &f);

// perpendicular Laplacian operator
const Field2D Delp2(const Field2D &f);
const Field3D Delp2(const Field3D &f, BoutReal zsmooth=-1.0);
const FieldPerp Delp2(const FieldPerp &f, BoutReal zsmooth=-1.0);

// Perpendicular Laplacian, keeping y derivatives
const Field2D Laplace_perp(const Field2D &f);
const Field3D Laplace_perp(const Field3D &f);

// Parallel Laplacian operator
const Field2D Laplace_par(const Field2D &f);
const Field3D Laplace_par(const Field3D &f);

// Full Laplacian operator (par + perp)
const Field2D Laplace(const Field2D &f);
const Field3D Laplace(const Field3D &f);

// Terms of form b0 x Grad(phi) dot Grad(A)
const Field2D b0xGrad_dot_Grad(const Field2D &phi, const Field2D &A);
const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field2D &A, CELL_LOC outloc=CELL_DEFAULT);
const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A);
const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc=CELL_DEFAULT);

// Poisson bracket methods
enum BRACKET_METHOD {BRACKET_STD=0, BRACKET_SIMPLE=1, BRACKET_ARAKAWA=2, BRACKET_CTU=3, BRACKET_ARAKAWA_OLD=4};
const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);
const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);
const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);
const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);

#endif /* __DIFOPS_H__ */
