/*!******************************************************************************
 * \file difops.hxx
 * 
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

/*!
 * Parallel derivative (central differencing) in Y
 * along unperturbed field 
 *
 * @param[in] var  The field to be differentiated
 * @param[in] outloc  The cell location where the output is needed (if staggered grids is enabled)
 * @param[in] method  The method to use. The default is set in the options. 
 */
const Field2D Grad_par(const Field2D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field2D Grad_par(const Field2D &var, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);
const Field3D Grad_par(const Field3D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field3D Grad_par(const Field3D &var, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);

/*!
 * Derivative along perturbed magnetic field
 * in Clebsch coordinate system
 *
 * b0 dot Grad  -  (1/B)b0 x Grad(apar) dot Grad
 *
 *
 * Combines the parallel and perpendicular calculation to include
 * grid-points at the corners.
 */
const Field3D Grad_parP(const Field3D &apar, const Field3D &f);

/*!
 * vpar times parallel derivative along unperturbed B-field (upwinding)
 *
 * \f[ 
 *    v\mathbf{b}_0 \cdot \nabla f
 * \f]
 * 
 * @param[in] v  The velocity in y direction
 * @param[in] f  The scalar field to be differentiated
 * 
 */
const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f);

/*!
 * vpar times parallel derivative along unperturbed B-field (upwinding)
 *
 * \f[ 
 *    v\mathbf{b}_0 \cdot \nabla f
 * \f]
 * 
 * @param[in] v  The velocity in y direction
 * @param[in] f  The scalar field to be differentiated
 * @param[in] outloc  The cell location of the output. By default this is the same as \p f
 * @param[in] method  The numerical method to use. The default is set in the options
 * 
 */
const Field3D Vpar_Grad_par(const Field &v, const Field &f, 
			    CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);

/*!
 * vpar times parallel derivative along unperturbed B-field (upwinding)
 *
 * \f[ 
 *    v\mathbf{b}_0 \cdot \nabla f
 * \f]
 * 
 * @param[in] v  The velocity in y direction
 * @param[in] f  The scalar field to be differentiated
 * @param[in] outloc  The cell location of the output. By default this is the same as \p f
 * @param[in] method  The numerical method to use. The default is set in the options
 * 
 */
const Field3D Vpar_Grad_par(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);


/*!
 * parallel divergence operator
 * 
 * \f[
 *  B \partial_{||}(f/B) = B \nabla\cdot (\mathbf{b}f/B )
 * \f]
 * 
 * @param[in] f  The component of a vector along the magnetic field 
 * 
 */
const Field2D Div_par(const Field2D &f);

/*!
 * parallel divergence operator
 * 
 * \f[
 *  B \partial_{||}(f/B) = B \nabla\cdot (\mathbf{b}f/B )
 * \f]
 * 
 * @param[in] f  The component of a vector along the magnetic field 
 * @param[in] outloc  The cell location for the result. By default the same as \p f
 * @param[in] method  The numerical method to use
 * 
 */
const Field3D Div_par(const Field3D &f, 
		      CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);

/*!
 * parallel divergence operator
 * 
 * \f[
 *  B \partial_{||}(f/B) = B \nabla\cdot (\mathbf{b}f/B )
 * \f]
 * 
 * @param[in] f  The component of a vector along the magnetic field 
 * @param[in] outloc  The cell location for the result. By default the same as \p f
 * @param[in] method  The numerical method to use
 * 
 */
const Field3D Div_par(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT);

// Flux methods. Model divergence of flux: df/dt =  Div(v * f)
const Field3D Div_par_flux(const Field3D &v, const Field3D &f, 
		      CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
const Field3D Div_par_flux(const Field3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT);

// Divergence of a parallel flow: Div(f*v)
// Both f and v are interpolated onto cell boundaries
// using 2nd order central difference, then multiplied together
// to get the flux at the boundary.
const Field3D Div_par(const Field3D &f, const Field3D &v);

/*!
 * second parallel derivative
 * \f[
 *    (\mathbf{b} dot \nabla)(\mathbf{b} dot \nabla)
 * \f]
 *
 * Note: For parallel Laplacian use LaplacePar
 */
const Field2D Grad2_par2(const Field2D &f);

/*!
 * second parallel derivative
 * \f[
 *    (\mathbf{b} dot \nabla)(\mathbf{b} dot \nabla)
 * \f]
 *
 * Note: For parallel Laplacian use LaplacePar
 *
 * @param[in] f The field to be differentiated
 * @param[in] outloc The cell location of the result 
 */
const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);

/*!
 * Parallel derivatives, converting between cell-centred and lower cell boundary
 * These are a simple way to do staggered differencing
 */
const Field3D Grad_par_CtoL(const Field3D &var);
const Field3D Vpar_Grad_par_LCtoC(const Field &v, const Field &f);
const Field3D Grad_par_LtoC(const Field3D &var);
const Field3D Div_par_LtoC(const Field2D &var);
const Field3D Div_par_LtoC(const Field3D &var);
const Field3D Div_par_CtoL(const Field2D &var);
const Field3D Div_par_CtoL(const Field3D &var);

/*!
 * Parallel divergence of diffusive flux, K*Grad_par
 * 
 * \f[
 *    \nabla \cdot ( \mathbf{b}_0 kY (\mathbf{b}_0 \cdot \nabla) f )
 * \f]
 *
 * @param[in] kY  The diffusion coefficient
 * @param[in] f   The field whose gradient drives a flux
 */
const Field2D Div_par_K_Grad_par(BoutReal kY, const Field2D &f);
const Field3D Div_par_K_Grad_par(BoutReal kY, const Field3D &f);
const Field2D Div_par_K_Grad_par(const Field2D &kY, const Field2D &f);
const Field3D Div_par_K_Grad_par(const Field2D &kY, const Field3D &f);
const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field2D &f);
const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field3D &f);

/*!
 * Perpendicular Laplacian operator
 *
 * This version only includes terms in X and Z, dropping
 * derivatives in Y. This is the inverse operation to 
 * the Laplacian inversion class. 
 *
 * For the full perpendicular Laplacian, use Laplace_perp
 */
const Field2D Delp2(const Field2D &f);
const Field3D Delp2(const Field3D &f, BoutReal zsmooth=-1.0);
const FieldPerp Delp2(const FieldPerp &f, BoutReal zsmooth=-1.0);

/*!
 * Perpendicular Laplacian, keeping y derivatives
 *
 * 
 */
const Field2D Laplace_perp(const Field2D &f);
const Field3D Laplace_perp(const Field3D &f);

/*!
 * Parallel Laplacian operator
 *
 */
const Field2D Laplace_par(const Field2D &f);
const Field3D Laplace_par(const Field3D &f);

/*!
 * Full Laplacian operator (par + perp)
 */
const Field2D Laplace(const Field2D &f);
const Field3D Laplace(const Field3D &f);

/*!
 * Terms of form b0 x Grad(phi) dot Grad(A)
 * 
 */
const Field2D b0xGrad_dot_Grad(const Field2D &phi, const Field2D &A);

/*!
 * Terms of form 
 *
 * \f[
 *   \mathbf{b}_0 \times \nabla \phi \cdot \nabla A
 * \f]
 * 
 * @param[in] phi The scalar potential
 * @param[in] A   The field being advected
 * @param[in] outloc  The cell location where the result is defined. By default the same as A.
 */
const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field2D &A, CELL_LOC outloc=CELL_DEFAULT);
const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A);
const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc=CELL_DEFAULT);

/*!
 * Poisson bracket methods
 */
enum BRACKET_METHOD {BRACKET_STD=0,   ///< Use b0xGrad_dot_Grad
                     BRACKET_SIMPLE=1, ///< Keep only terms in X-Z
                     BRACKET_ARAKAWA=2, ///< Arakawa method in X-Z (optimised)
                     BRACKET_CTU=3, ///< Corner Transport Upwind (CTU) method. Explicit method only, needs the timestep from the solver
                     BRACKET_ARAKAWA_OLD=4, ///< Older version, for regression testing of optimised version.
                     BRACKET_ARAKAWA_SDI=5 ///< SingleDataIterator version for testing
};

/*!
 * Compute advection operator terms, which can be cast as
 * antisymmetric Poisson brackets
 *
 * \f[
 *   [f, g] = (1/B) \mathbf{b}_0 \times \nabla f  \cdot \nabla g
 * \f]
 * 
 * @param[in] f  The potential
 * @param[in] f  The field being advected
 * @param[in] method   The method to use
 * @param[in] outloc   The cell location where the result is defined. Default is the same as g
 * @param[in] solver   Pointer to the time integration solver
 * 
 */
const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);
const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);
const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);
const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc=CELL_DEFAULT, Solver *solver = NULL);

#endif /* __DIFOPS_H__ */
