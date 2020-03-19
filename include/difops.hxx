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

#include "bout/deprecated.hxx"
#include "bout/solver.hxx"

/*!
 * Parallel derivative (central differencing) in Y
 * along unperturbed field
 *
 * @param[in] var  The field to be differentiated
 * @param[in] outloc  The cell location where the output is needed (if staggered grids is
 * enabled)
 * @param[in] method  The method to use. The default is set in the options.
 */
const Coordinates::metric_field_type Grad_par(const Field2D& var,
                                              CELL_LOC outloc = CELL_DEFAULT,
                                              const std::string& method = "DEFAULT");
DEPRECATED(const Coordinates::metric_field_type Grad_par(const Field2D& var,
                                                         const std::string& method,
                                                         CELL_LOC outloc = CELL_DEFAULT));
inline const Coordinates::metric_field_type Grad_par(const Field2D& var, CELL_LOC outloc,
                                                     DIFF_METHOD method) {
  return Grad_par(var, outloc, toString(method));
}
DEPRECATED(inline const Coordinates::metric_field_type Grad_par(const Field2D& var,
                                                                DIFF_METHOD method,
                                                                CELL_LOC outloc)) {
  return Grad_par(var, outloc, toString(method));
}

const Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
                       const std::string& method = "DEFAULT");
DEPRECATED(const Field3D Grad_par(const Field3D& var, const std::string& method,
                                  CELL_LOC outloc = CELL_DEFAULT));
inline const DEPRECATED(Field3D Grad_par(const Field3D& var, CELL_LOC outloc,
                                         DIFF_METHOD method)) {
  return Grad_par(var, outloc, toString(method));
};
DEPRECATED(inline const
    Field3D Grad_par(const Field3D& var, DIFF_METHOD method, CELL_LOC outloc)) {
  return Grad_par(var, outloc, toString(method));
}

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
const Field3D Grad_parP(const Field3D& apar, const Field3D& f);

/*!
 * vpar times parallel derivative along unperturbed B-field (upwinding)
 *
 * \f[
 *    v\mathbf{b}_0 \cdot \nabla f
 * \f]
 * 
 *
 * @param[in] v  The velocity in y direction
 * @param[in] f  The scalar field to be differentiated
 * @param[in] outloc  The cell location of the output. By default this is the same as \p f
 * @param[in] method  The numerical method to use. The default is set in the options
 *
 */
const Coordinates::metric_field_type Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                                   CELL_LOC outloc = CELL_DEFAULT,
                                                   const std::string& method = "DEFAULT");
DEPRECATED(const Coordinates::metric_field_type Vpar_Grad_par(
    const Field2D& v, const Field2D& f, const std::string& method,
    CELL_LOC outloc = CELL_DEFAULT));
inline const Coordinates::metric_field_type
Vpar_Grad_par(const Field2D& v, const Field2D& f, CELL_LOC outloc, DIFF_METHOD method) {
  return Vpar_Grad_par(v, f, outloc, toString(method));
}
DEPRECATED(inline const Coordinates::metric_field_type Vpar_Grad_par(const Field2D& v,
                                                                     const Field2D& f,
                                                                     DIFF_METHOD method,
                                                                     CELL_LOC outloc)) {
  return Vpar_Grad_par(v, f, outloc, toString(method));
}

const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                            CELL_LOC outloc = CELL_DEFAULT,
                            const std::string& method = "DEFAULT");
DEPRECATED(const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                                       const std::string& method,
                                       CELL_LOC outloc = CELL_DEFAULT));
inline const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
                                   DIFF_METHOD method) {
  return Vpar_Grad_par(v, f, outloc, toString(method));
}
DEPRECATED(inline const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                                              DIFF_METHOD method, CELL_LOC outloc)) {
  return Vpar_Grad_par(v, f, outloc, toString(method));
}

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
const Coordinates::metric_field_type Div_par(const Field2D& f,
                                             CELL_LOC outloc = CELL_DEFAULT,
                                             const std::string& method = "DEFAULT");
DEPRECATED(const Coordinates::metric_field_type Div_par(const Field2D& f,
                                                        const std::string& method,
                                                        CELL_LOC outloc = CELL_DEFAULT));
inline const Coordinates::metric_field_type Div_par(const Field2D& f, CELL_LOC outloc,
                                                    DIFF_METHOD method) {
  return Div_par(f, outloc, toString(method));
}
DEPRECATED(inline const Coordinates::metric_field_type Div_par(const Field2D& f,
                                                               DIFF_METHOD method,
                                                               CELL_LOC outloc)) {
  return Div_par(f, outloc, toString(method));
}

const Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& method = "DEFAULT");
DEPRECATED(const Field3D Div_par(const Field3D& f, const std::string& method,
                                 CELL_LOC outloc = CELL_DEFAULT));
inline const Field3D Div_par(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
  return Div_par(f, outloc, toString(method));
}
DEPRECATED(inline const Field3D Div_par(const Field3D& f, DIFF_METHOD method,
                                        CELL_LOC outloc)) {
  return Div_par(f, outloc, toString(method));
}

// Divergence of a parallel flow: Div(f*v)
// Both f and v are interpolated onto cell boundaries
// using 2nd order central difference, then multiplied together
// to get the flux at the boundary.
const Field3D Div_par(const Field3D& f, const Field3D& v);

// Flux methods. Model divergence of flux: df/dt =  Div(v * f)
// TODO : Should we add Field2D versions?
const Field3D Div_par_flux(const Field3D& v, const Field3D& f,
                           CELL_LOC outloc = CELL_DEFAULT,
                           const std::string& method = "DEFAULT");
DEPRECATED(const Field3D Div_par_flux(const Field3D& v, const Field3D& f,
                                      const std::string& method,
                                      CELL_LOC outloc = CELL_DEFAULT));
inline const Field3D Div_par_flux(const Field3D& v, const Field3D& f, CELL_LOC outloc,
                                  DIFF_METHOD method) {
  return Div_par_flux(v, f, outloc, toString(method));
}
DEPRECATED(inline const Field3D Div_par_flux(const Field3D& v, const Field3D& f,
                                             DIFF_METHOD method,
                                             CELL_LOC outloc = CELL_DEFAULT)) {
  return Div_par_flux(v, f, outloc, toString(method));
}

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
const Coordinates::metric_field_type Grad2_par2(const Field2D& f,
                                                CELL_LOC outloc = CELL_DEFAULT,
                                                const std::string& method = "DEFAULT");
inline const Coordinates::metric_field_type Grad2_par2(const Field2D& f, CELL_LOC outloc,
                                                       DIFF_METHOD method) {
  return Grad2_par2(f, outloc, toString(method));
}

const Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");
inline const Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
  return Grad2_par2(f, outloc, toString(method));
}

/*!
 * Parallel derivatives, converting between cell-centred and lower cell boundary
 * These are a simple way to do staggered differencing
 */
[[deprecated(
    "Grad_par_CtoL is deprecated. Staggering is now supported in Grad_par.")]]
inline const Field3D Grad_par_CtoL(const Field3D &var) {
  ASSERT2(var.getLocation() == CELL_CENTRE);
  return Grad_par(var, CELL_YLOW);
}
[[deprecated(
    "Grad_par_CtoL is deprecated. Staggering is now supported in Grad_par.")]]
inline const Coordinates::metric_field_type Grad_par_CtoL(const Field2D &var) {
  ASSERT2(var.getLocation() == CELL_CENTRE);
  return Grad_par(var, CELL_YLOW);
}
[[deprecated(
    "Vpar_Grad_par_LCtoC is deprecated. Staggering is now supported in Vpar_Grad_par.")]]
inline const Field3D Vpar_Grad_par_LCtoC(const Field3D& v, const Field3D& f,
    const std::string& region="RGN_NOBNDRY") {
  ASSERT2(v.getLocation() == CELL_YLOW);
  ASSERT2(f.getLocation() == CELL_CENTRE);
  return Vpar_Grad_par(v, f, CELL_CENTRE, region);
}
[[deprecated(
    "Vpar_Grad_par_LCtoC is deprecated. Staggering is now supported in Vpar_Grad_par.")]]
inline const Field3D Vpar_Grad_par_LCtoC(const Field3D& v, const Field3D& f,
    REGION region=RGN_NOBNDRY) {
  ASSERT2(v.getLocation() == CELL_YLOW);
  ASSERT2(f.getLocation() == CELL_CENTRE);
  return Vpar_Grad_par(v, f, CELL_CENTRE, toString(region));
}
[[deprecated(
    "Grad_par_LtoC is deprecated. Staggering is now supported in Grad_par.")]]
inline const Field3D Grad_par_LtoC(const Field3D &var) {
  ASSERT2(var.getLocation() == CELL_YLOW);
  return Grad_par(var, CELL_CENTRE);
}
[[deprecated(
    "Grad_par_LtoC is deprecated. Staggering is now supported in Grad_par.")]]
inline const Coordinates::metric_field_type Grad_par_LtoC(const Field2D &var) {
  ASSERT2(var.getLocation() == CELL_YLOW);
  return Grad_par(var, CELL_CENTRE);
}
[[deprecated(
    "Div_par_LtoC is deprecated. Staggering is now supported in Grad_par.")]]
inline const Field3D Div_par_LtoC(const Field3D &var) {
  ASSERT2(var.getLocation() == CELL_YLOW);
  return Div_par(var, CELL_CENTRE);
}
[[deprecated(
    "Div_par_LtoC is deprecated. Staggering is now supported in Grad_par.")]]
inline const Coordinates::metric_field_type Div_par_LtoC(const Field2D &var) {
  ASSERT2(var.getLocation() == CELL_YLOW);
  return Div_par(var, CELL_CENTRE);
}
[[deprecated(
    "Div_par_CtoL is deprecated. Staggering is now supported in Grad_par.")]]
inline const Field3D Div_par_CtoL(const Field3D &var) {
  ASSERT2(var.getLocation() == CELL_CENTRE);
  return Div_par(var, CELL_YLOW);
}
[[deprecated(
    "Div_par_CtoL is deprecated. Staggering is now supported in Grad_par.")]]
inline const Coordinates::metric_field_type Div_par_CtoL(const Field2D &var) {
  ASSERT2(var.getLocation() == CELL_CENTRE);
  return Div_par(var, CELL_YLOW);
}

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
const Coordinates::metric_field_type Div_par_K_Grad_par(BoutReal kY, const Field2D& f,
                                                        CELL_LOC outloc = CELL_DEFAULT);
const Field3D Div_par_K_Grad_par(BoutReal kY, const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
const Coordinates::metric_field_type
Div_par_K_Grad_par(const Field2D& kY, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
const Field3D Div_par_K_Grad_par(const Field2D &kY, const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field2D &f, CELL_LOC outloc=CELL_DEFAULT);
const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);

/*!
 * Perpendicular Laplacian operator
 *
 * This version only includes terms in X and Z, dropping
 * derivatives in Y. This is the inverse operation to 
 * the Laplacian inversion class. 
 *
 * For the full perpendicular Laplacian, use Laplace_perp
 */
const Coordinates::metric_field_type Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
const Field3D Delp2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
const FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT,
                      bool useFFT = true);

/*!
 * Perpendicular Laplacian, keeping y derivatives
 *
 * 
 */
const Coordinates::metric_field_type Laplace_perp(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                                                  const std::string& dfdy_boundary_condition = "free_o3",
                                                  const std::string& dfdy_region = "");
const Field3D Laplace_perp(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                           const std::string& dfdy_boundary_condition = "free_o3",
                           const std::string& dfdy_region = "");

/*!
 * Parallel Laplacian operator
 *
 */
const Coordinates::metric_field_type Laplace_par(const Field2D& f,
                                                 CELL_LOC outloc = CELL_DEFAULT);
const Field3D Laplace_par(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);

/*!
 * Full Laplacian operator (par + perp)
 */
const Coordinates::metric_field_type Laplace(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                                             const std::string& dfdy_boundary_condition = "free_o3",
                                             const std::string& dfdy_region = "");
const Field3D Laplace(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& dfdy_boundary_condition = "free_o3",
                      const std::string& dfdy_region = "");

/*!
 * Inverse of Laplacian operator in LaplaceXY solver
 */
const Field2D Laplace_perpXY(const Field2D& A, const Field2D& f);

/*!
 * Terms of form b0 x Grad(phi) dot Grad(A)
 * 
 */
const Coordinates::metric_field_type
b0xGrad_dot_Grad(const Field2D& phi, const Field2D& A, CELL_LOC outloc = CELL_DEFAULT);

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
const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A, CELL_LOC outloc=CELL_DEFAULT);
const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc=CELL_DEFAULT);


/*!
 * X-Z Finite Volume diffusion operator
 */
const Field3D Div_Perp_Lap_FV(const Field3D &a, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);


/*!
 * Poisson bracket methods
 */
enum class BRACKET_METHOD {
  standard,   ///< Use b0xGrad_dot_Grad
  simple,     ///< Keep only terms in X-Z
  arakawa,    ///< Arakawa method in X-Z (optimised)
  ctu,        ///< Corner Transport Upwind (CTU) method. Explicit method only, needs the
              ///  timestep from the solver
  arakawa_old ///< Older version, for regression testing of optimised version.
};
constexpr BRACKET_METHOD BRACKET_STD = BRACKET_METHOD::standard;
constexpr BRACKET_METHOD BRACKET_SIMPLE = BRACKET_METHOD::simple;
constexpr BRACKET_METHOD BRACKET_ARAKAWA = BRACKET_METHOD::arakawa;
constexpr BRACKET_METHOD BRACKET_CTU = BRACKET_METHOD::ctu;
constexpr BRACKET_METHOD BRACKET_ARAKAWA_OLD = BRACKET_METHOD::arakawa_old;

/*!
 * Compute advection operator terms, which can be cast as
 * antisymmetric Poisson brackets
 *
 * \f[
 *   [f, g] = (1/B) \mathbf{b}_0 \times \nabla f  \cdot \nabla g
 * \f]
 * 
 * @param[in] f  The potential
 * @param[in] g  The field being advected
 * @param[in] method   The method to use
 * @param[in] outloc   The cell location where the result is defined. Default is the same as g
 * @param[in] solver   Pointer to the time integration solver
 * 
 */
const Coordinates::metric_field_type bracket(const Field2D& f, const Field2D& g,
                                             BRACKET_METHOD method = BRACKET_STD,
                                             CELL_LOC outloc = CELL_DEFAULT,
                                             Solver* solver = nullptr);
const Field3D bracket(const Field2D &f, const Field3D &g,
                      BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc = CELL_DEFAULT,
                      Solver *solver = nullptr);
const Field3D bracket(const Field3D &f, const Field2D &g,
                      BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc = CELL_DEFAULT,
                      Solver *solver = nullptr);
const Field3D bracket(const Field3D &f, const Field3D &g,
                      BRACKET_METHOD method = BRACKET_STD, CELL_LOC outloc = CELL_DEFAULT,
                      Solver *solver = nullptr);

#endif /* __DIFOPS_H__ */
