/*
  Finite-volume discretisation methods. Flux-conservative form
 */

#ifndef BOUT_FV_OPS_H
#define BOUT_FV_OPS_H

#include "bout/bout_types.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/output_bout_types.hxx" // NOLINT(unused-includes, misc-include-cleaner)
#include "bout/vector2d.hxx"

namespace FV {
/// Vorticity:
/// \f[
///   \nabla (a \cdot \nabla_\perp(f))
/// \f]
// Div ( a Grad_perp(f) ) -- ∇⊥ ( a ⋅ ∇⊥ f) -- Vorticity
Field3D Div_a_Grad_perp(const Field3D& a, const Field3D& f);

[[deprecated("Please use Div_a_Grad_perp instead")]] inline Field3D
Div_a_Laplace_perp(const Field3D& a, const Field3D& f) {
  return Div_a_Grad_perp(a, f);
}

/// Divergence of a parallel diffusion
/// \f[
///   \nabla_\parallel(k \cdot \nabla_\parallel(f) )
/// \f]
Field3D Div_par_K_Grad_par(const Field3D& k, const Field3D& f, bool bndry_flux = true);

/// 4th-order derivative in Y, using derivatives
/// on cell boundaries.
///
/// A one-sided 3rd-order derivative, given a value
/// at a boundary is:
///
/// \f[
/// \frac{d^3f}{dx^3} \simeq \tfrac{16}{5} f_b - 6 f_0 + 4 f_1 - \tfrac{6}{5} f_2
/// \f]
///
/// where:
///
/// - \f$f_b\f$ is the value on the boundary,
/// - \f$f_0\f$ is the cell to the left of the boundary,
/// - \f$f_1\f$ to the left of \f$f_0\f$, and
/// - \f$f_2\f$ to the left of f_1:
///
/// .. code:: text
///    f_2 | f_1 | f_0 |
///                   f_b
///
/// NB: Uses to/from FieldAligned coordinates
///
/// No fluxes through domain boundaries

Field3D D4DY4(const Field3D& d, const Field3D& f);

/// 4th-order dissipation term
///
/// \f[
/// \frac{d^3f}{dx^3} \simeq \tfrac{16}{5} f_b - 6 f_0 + 4 f_1 - \tfrac{6}{5} f_2
/// \f]
///
/// where:
///
/// - \f$f_b\f$ is the value on the boundary,
/// - \f$f_0\f$ is the cell to the left of the boundary,
/// - \f$f_1\f$ to the left of \f$f_0\f$, and
/// - \f$f_2\f$ to the left of f_1:
///
/// .. code:: text
///    f_2 | f_1 | f_0 |
///                   f_b
Field3D D4DY4_Index(const Field3D& f, bool bndry_flux = true);

// Forward declarations of flux limiters
// If you want to use your own flux limiter, you need to
// #include <bout/fv_ops_impl.hxx> to instantiate the templates.
struct Upwind;
struct Fromm;
struct MinMod;
struct MC;
struct Superbee;
struct VanAlbada;
struct WENO3;

/*!
   * Communicate fluxes between processors
   * Takes values in guard cells, and adds them to cells
   */
void communicateFluxes(Field3D& f);

/// Finite volume parallel divergence
///
/// Preserves the sum of f*J*dx*dy*dz over the domain
///
/// @param[in] f_in   The field being advected.
///                   This will be reconstructed at cell faces
///                   using the given CellEdges method
/// @param[in] v_in   The advection velocity.
///                   This will be interpolated to cell boundaries
///                   using linear interpolation
/// @param[in] wave_speed_in  Local maximum speed of all waves in the system at each
//                            point in space
/// @param[in] fixflux     Fix the flux at the boundary to be the value at the
///                        midpoint (for boundary conditions)
///
/// NB: Uses to/from FieldAligned coordinates
template <typename CellEdges = MC>
Field3D Div_par(const Field3D& f_in, const Field3D& v_in, const Field3D& wave_speed_in,
                bool fixflux = true);

/*!
   * Div ( n * v )  -- Magnetic drifts
   *
   * This uses the expression
   *
   * Div( A ) = 1/J * d/di ( J * A^i )
   *
   * Hence the input vector should be contravariant
   *
   * Note: Uses to/from FieldAligned
   *
   */
template <typename CellEdges = MC>
Field3D Div_f_v(const Field3D& n_in, const Vector3D& v, bool bndry_flux);

/*!
   * X-Z Finite Volume diffusion operator
   */
Field3D Div_Perp_Lap(const Field3D& a, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT);

/// Finite volume parallel divergence
///
/// NOTE: Modified version, applies limiter to velocity and field
///       Performs better (smaller overshoots) than Div_par
///
/// Preserves the sum of f*J*dx*dy*dz over the domain
///
/// @param[in] f_in   The field being advected.
///                   This will be reconstructed at cell faces
///                   using the given CellEdges method
/// @param[in] v_in   The advection velocity.
///                   This will be interpolated to cell boundaries
///                   using linear interpolation
/// @param[in] wave_speed_in  Local maximum speed of all waves in the system at each
//                            point in space
/// @param[in] fixflux     Fix the flux at the boundary to be the value at the
///                        midpoint (for boundary conditions)
///
/// @param[out] flow_ylow    Flow at the lower Y cell boundary
///                          Already includes area factor * flux
///                          For FCI fields this diagnostic is currently set to zero.
template <typename CellEdges = MC>
Field3D Div_par_mod(const Field3D& f_in, const Field3D& v_in,
                    const Field3D& wave_speed_in, Field3D& flow_ylow, bool fixflux = true,
                    bool dissipative = false);

/// This operator calculates Div_par(f v v)
/// It is used primarily (only?) in the parallel momentum equation.
///
/// This operator is used rather than Div(f fv) so that the values of
/// f and v are consistent with other advection equations: The product
/// fv is not interpolated to cell boundaries.
template <typename CellEdges = MC>
Field3D Div_par_fvv(const Field3D& f_in, const Field3D& v_in,
                    const Field3D& wave_speed_in, bool fixflux = true);

/// Calculates viscous heating due to numerical momentum fluxes
/// and flow of kinetic energy (in flow_ylow)
template <typename CellEdges = MC>
Field3D Div_par_fvv_heating(const Field3D& f_in, const Field3D& v_in,
                            const Field3D& wave_speed_in, Field3D& flow_ylow,
                            bool fixflux = true);

/// Div ( a g Grad_perp(f) )  -- Perpendicular gradient-driven advection
///
/// This version uses a slope limiter to calculate cell edge values of g in X,
/// the advects the upwind cell edge.
///
/// 1st order upwinding is used in Y.
template <typename CellEdges = MC>
Field3D Div_a_Grad_perp_limit(const Field3D& a, const Field3D& g, const Field3D& f);

/// Div ( a g Grad_perp(f) )
///
/// This version uses pre-computed coefficient. It can also be used
class dagp_fv {
public:
  Field3D operator()(const Field3D& a, const Field3D& f, Field3D& low_xlow,
                     Field3D& flow_zlow, bool upwinding);
  Field3D operator()(const Field3D& a, const Field3D& f, bool upwinding);
  dagp_fv(Mesh& mesh);
  dagp_fv& operator*=(BoutReal fac);
  dagp_fv& operator/=(BoutReal fac);

private:
  template <bool extra, bool upwinding>
  Field3D operator()(const Field3D& a, const Field3D& f, Field3D* low_xlow,
                     Field3D* flow_zlow);
  Field3D fac_XX;
  Field3D fac_XZ;
  Field3D fac_ZX;
  Field3D fac_ZZ;
  Field3D volume;

  bool isNormalised{false};

  template <bool upwinding>
  BoutReal xflux(const Field3D& a, const Field3D& f, const Ind3D& i);
  template <bool upwinding>
  BoutReal zflux(const Field3D& a, const Field3D& f, const Ind3D& i);
};

std::shared_ptr<dagp_fv> getDagp_fv(Mesh* mesh, BoutReal rho_s0);

} // namespace FV
#endif // BOUT_FV_OPS_H
