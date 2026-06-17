/**************************************************************************
 * Shared helpers for SUNDIALS N_Vector backends used by BOUT++ solvers.
 *
 * Copyright 2010-2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_SUNDIALS_NVECTOR_INTERFACE_HXX
#define BOUT_SUNDIALS_NVECTOR_INTERFACE_HXX

#include "bout/bout_enum_class.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/solver.hxx"
#include "bout/sundials_backports.hxx"

#if BOUT_HAS_SUNDIALS_MANYVECTOR
#include "nvector.hxx"
#endif

#include <algorithm>
#include <iterator>
#include <vector>

/// Supported SUNDIALS ``N_Vector`` backends that can be selected at runtime.
BOUT_ENUM_CLASS(NVectorType, Sundials, ManyVector);

/// Helper for moving solver state between BOUT++ fields and a chosen ``N_Vector``
/// backend without exposing backend-specific details in each solver wrapper.
class SundialsNVectorInterface {
public:
  SundialsNVectorInterface(Solver& solver_in, const sundials::Context& ctx_in,
                           NVectorType nvector_type_in)
      : solver(solver_in), ctx(ctx_in), nvector_type(nvector_type_in) {}

  /// Return ``true`` when using the field-backed ManyVector implementation.
  bool use_manyvector() const { return nvector_type == NVectorType::ManyVector; }

  /// Human-readable backend name for logs and error messages.
  const char* backend_name() const {
    return use_manyvector() ? "SUNDIALS ManyVector-backed custom" : "SUNDIALS parallel";
  }

  /// Throw if the selected backend requires SUNDIALS ManyVector support that was
  /// not enabled when BOUT++ was built.
  void ensure_manyvector_available() const {
    if (use_manyvector() and not BOUT_HAS_SUNDIALS_MANYVECTOR) {
      throw BoutException("SUNDIALS ManyVector backend requested, but BOUT++ was "
                          "built without SUNDIALS ManyVector support");
    }
  }

  /// Create the main solver state vector using the selected backend.
  N_Vector create_state_vector(int local_N, int neq) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      auto* vec = nvector_from_state();
      if (vec == nullptr) {
        throw BoutException("BOUT N_Vector failed\n");
      }
      return vec;
#endif
    }

    auto* vec = callWithSUNContext(N_VNew_Parallel, ctx, BoutComm::get(), local_N, neq);
    if (vec == nullptr) {
      throw BoutException("SUNDIALS memory allocation failed\n");
    }
    solver.save_vars(N_VGetArrayPointer(vec));
    return vec;
  }

  /// Clone a vector of the same backend/layout, with a descriptive error on failure.
  N_Vector clone_vector_like(N_Vector source, const char* description) const {
    auto* vec = N_VClone(source);
    if (vec == nullptr) {
      throw BoutException("{}\n", description);
    }
    return vec;
  }

  /// Copy the current solution state from ``u`` into the solver fields.
  void copy_state_from_vector(N_Vector u) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_state(u);
      return;
#endif
    }

    solver.load_vars(N_VGetArrayPointer(u));
  }

  /// Copy the solver fields into the solution vector ``u``.
  void copy_state_to_vector(N_Vector u) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_state(u);
      return;
#endif
    }

    solver.save_vars(N_VGetArrayPointer(u));
  }

  /// Copy time-derivative values from ``du`` into the solver derivative fields.
  void copy_deriv_from_vector(N_Vector du) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_deriv(du);
      return;
#endif
    }

    solver.load_derivs(N_VGetArrayPointer(du));
  }

  /// Copy the solver derivative fields into ``du``.
  void copy_deriv_to_vector(N_Vector du) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_deriv(du);
      return;
#endif
    }

    solver.save_derivs(N_VGetArrayPointer(du));
  }

  /// Fill every evolved entry in ``vec`` with per-field constants, respecting
  /// whether each field evolves boundary cells.
  void fill_vector_values(N_Vector vec, const std::vector<BoutReal>& f2d_values,
                          const std::vector<BoutReal>& f3d_values) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      std::size_t i = 0;
      for (std::size_t j = 0; j < f2d_values.size(); ++j, ++i) {
        fill_field(BoutNVector::get<Field2D>(vec, i), f2d_values[j],
                   solver.f2d[j].evolve_bndry);
      }
      for (std::size_t j = 0; j < f3d_values.size(); ++j, ++i) {
        fill_field(BoutNVector::get<Field3D>(vec, i), f3d_values[j],
                   solver.f3d[j].evolve_bndry);
      }
      return;
#endif
    }

    BoutReal* option_data = N_VGetArrayPointer(vec);
    int p = 0;
    for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_BNDRY")) {
      fill_vector_values_op(i2d, option_data, p, f2d_values, f3d_values, true);
    }
    for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_NOBNDRY")) {
      fill_vector_values_op(i2d, option_data, p, f2d_values, f3d_values, false);
    }
  }

  /// Fill the IDA ``id`` vector with ``1`` for differential variables and ``0``
  /// for constrained algebraic variables.
  void fill_id_vector(N_Vector vec) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      std::size_t i = 0;
      for (const auto& field : solver.f2d) {
        fill_field(BoutNVector::get<Field2D>(vec, i), field.constraint ? 0.0 : 1.0,
                   field.evolve_bndry);
        ++i;
      }
      for (const auto& field : solver.f3d) {
        fill_field(BoutNVector::get<Field3D>(vec, i), field.constraint ? 0.0 : 1.0,
                   field.evolve_bndry);
        ++i;
      }
      return;
#endif
    }

    solver.set_id(N_VGetArrayPointer(vec));
  }

  /// Apply the IDA residual update ``residual -= du`` only on differential
  /// components selected by ``id``.
  void subtract_differential_term(N_Vector residual, N_Vector du, N_Vector id) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      std::size_t i = 0;
      for (const auto& field : solver.f2d) {
        subtract_field(BoutNVector::get<Field2D>(residual, i),
                       BoutNVector::get<Field2D>(du, i), BoutNVector::get<Field2D>(id, i),
                       field.evolve_bndry);
        ++i;
      }
      for (const auto& field : solver.f3d) {
        subtract_field(BoutNVector::get<Field3D>(residual, i),
                       BoutNVector::get<Field3D>(du, i), BoutNVector::get<Field3D>(id, i),
                       field.evolve_bndry);
        ++i;
      }
      return;
#endif
    }

    const auto length = N_VGetLocalLength_Parallel(id);
    const BoutReal* id_data = N_VGetArrayPointer(id);
    BoutReal* residual_data = N_VGetArrayPointer(residual);
    const BoutReal* du_data = N_VGetArrayPointer(du);
    for (int i = 0; i < length; ++i) {
      if (id_data[i] > 0.5) {
        residual_data[i] -= du_data[i];
      }
    }
  }

private:
  Solver& solver;
  const sundials::Context& ctx;
  NVectorType nvector_type;

  void fill_vector_values_op(Ind2D UNUSED(i2d), BoutReal* option_data, int& p,
                             const std::vector<BoutReal>& f2d_values,
                             const std::vector<BoutReal>& f3d_values, bool bndry) const {
    for (std::size_t i = 0; i < f2d_values.size(); ++i) {
      if (bndry && !solver.f2d[i].evolve_bndry) {
        continue;
      }
      option_data[p] = f2d_values[i];
      ++p;
    }

    for (int jz = 0; jz < bout::globals::mesh->LocalNz; ++jz) {
      for (std::size_t i = 0; i < f3d_values.size(); ++i) {
        if (bndry && !solver.f3d[i].evolve_bndry) {
          continue;
        }
        option_data[p] = f3d_values[i];
        ++p;
      }
    }
  }

#if BOUT_HAS_SUNDIALS_MANYVECTOR
  /// Build a field-backed ManyVector that aliases the solver state variables.
  N_Vector nvector_from_state() const {
    std::vector<N_Vector> subvectors;
    subvectors.reserve(solver.f2d.size() + solver.f3d.size());
    const auto inserter = std::back_inserter(subvectors);

    const auto field_to_nvector = [this](const auto& var_str) {
      return BoutNVector::create(ctx, *var_str.var, var_str.evolve_bndry);
    };
    std::transform(solver.f2d.cbegin(), solver.f2d.cend(), inserter, field_to_nvector);
    std::transform(solver.f3d.cbegin(), solver.f3d.cend(), inserter, field_to_nvector);
    return BoutNVector::create(ctx, subvectors);
  }

  /// Swap the active solver state fields with the subvectors contained in ``u``.
  void swap_state(const N_Vector u) const {
    std::size_t i = 0;
    for (auto& var_str : solver.f2d) {
      BoutNVector::swap(u, *var_str.var, i);
      ++i;
    }
    for (auto& var_str : solver.f3d) {
      BoutNVector::swap(u, *var_str.var, i);
      ++i;
    }
  }

  /// Swap the active solver derivative fields with the subvectors in ``du``.
  void swap_deriv(const N_Vector du) const {
    std::size_t i = 0;
    for (auto& var_str : solver.f2d) {
      BoutNVector::swap(du, *var_str.F_var, i);
      ++i;
    }
    for (auto& var_str : solver.f3d) {
      BoutNVector::swap(du, *var_str.F_var, i);
      ++i;
    }
  }

  template <typename FieldType>
  static void fill_field(FieldType& field, BoutReal value, bool evolve_bndry) {
    BOUT_FOR(i, field.getRegion(evolve_bndry ? RGN_ALL : RGN_NOBNDRY)) {
      field[i] = value;
    }
  }

  template <typename FieldType>
  static void subtract_field(FieldType& residual, const FieldType& du,
                             const FieldType& id, bool evolve_bndry) {
    BOUT_FOR(i, residual.getRegion(evolve_bndry ? RGN_ALL : RGN_NOBNDRY)) {
      if (id[i] > 0.5) {
        residual[i] -= du[i];
      }
    }
  }
#endif
};

#endif
